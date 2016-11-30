import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.manifold.t_sne import _joint_probabilities, _joint_probabilities_nn
from scipy.spatial.distance import squareform
from sklearn.neighbors import NearestNeighbors

def magic(data, kernel='gaussian', n_pca_components=None, t=8, knn=20, 
          knn_autotune=0, epsilon=0, rescale=True, k_knn=100, perplexity=30):

    if kernel not in ['gaussian', 'tsne']:
        raise RuntimeError('Invalid kernel type. Must be either "gaussian" or "tsne".')

    if n_pca_components != None:
        pca_loadings = run_pca(data, n_components=n_pca_components)
        pca_projected_data = np.dot(data, pca_loadings)
    else:
        pca_projected_data = data

    if kernel == 'gaussian':
        #run diffusion maps to get markov matrix
        L = compute_markov(pca_projected_data, knn=knn, epsilon=epsilon, 
                                       distance_metric='euclidean', knn_autotune=knn_autotune)

    else:
        #tsne kernel
        distances = pairwise_distances(pca_projected_data, squared=True)
        if k_knn > 0:
            neighbors_nn = np.argsort(distances, axis=0)[:, :k_knn]
            P = _joint_probabilities_nn(distances, neighbors_nn, perplexity, 1)
        else:
            P = _joint_probabilities(distances, perplexity, 1)
        P = squareform(P)

        #markov normalize P
        L = np.divide(P, np.sum(P, axis=1))

    #get imputed data matrix
    new_data, L_t = impute_fast(data, L, t, rescale_to_max=rescale)

    return new_data


def run_pca(data, n_components=100):

    X = data
    # Make sure data is zero mean
    X = np.subtract(X, np.amin(X))
    X = np.divide(X, np.amax(X))

    # Compute covariance matrix
    if (X.shape[1] < X.shape[0]):
        C = np.cov(X, rowvar=0)
    # if N>D, we better use this matrix for the eigendecomposition
    else:
        C = np.multiply((1/X.shape[0]), np.dot(X, X.T))

    # Perform eigendecomposition of C
    C[np.where(np.isnan(C))] = 0
    C[np.where(np.isinf(C))] = 0
    l, M = np.linalg.eig(C)

    # Sort eigenvectors in descending order
    ind = np.argsort(l)[::-1]
    l = l[ind]
    if n_components < 1:
        n_components = np.where(np.cumsum(np.divide(l, np.sum(l)), axis=0) >= n_components)[0][0] + 1
        print('Embedding into ' + str(n_components) + ' dimensions.')
    if n_components > M.shape[1]:
        n_components = M.shape[1]
        print('Target dimensionality reduced to ' + str(n_components) + '.')

    M = M[:, ind[:n_components]]
    l = l[:n_components]

    # Apply mapping on the data
    if X.shape[1] >= X.shape[0]:
        M = np.multiply(np.dot(X.T, M), (1 / np.sqrt(X.shape[0] * l)).T)

    return M


def impute_fast(data, L, t, rescale_to_max, L_t=None, tprev=None):

    #convert L to full matrix
    if issparse(L):
        L = L.todense()

    #L^t
    print('MAGIC: L_t = L^t')
    if L_t == None:
        L_t = np.linalg.matrix_power(L, t)
    else:
        L_t = np.dot(L_t, np.linalg.matrix_power(L, t-tprev))

    print('MAGIC: data_new = L_t * data')
    data_new = np.array(np.dot(L_t, data))

    #rescale data to 99th percentile
    if rescale_to_max == True:
        M99 = np.percentile(data, 99, axis=0)
        M100 = data.max(axis=0)
        indices = np.where(M99 == 0)[0]
        M99[indices] = M100[indices]
        M99_new = np.percentile(data_new, 99, axis=0)
        M100_new = data_new.max(axis=0)
        indices = np.where(M99_new == 0)[0]
        M99_new[indices] = M100_new[indices]
        max_ratio = np.divide(M99, M99_new)
        data_new = np.multiply(data_new, np.matlib.repmat(max_ratio, len(data), 1))
    
    return data_new, L_t


def compute_markov(data, knn=10, epsilon=1, distance_metric='euclidean', knn_autotune=0):

    N = data.shape[0]

    # Nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=knn, metric=distance_metric).fit(data)
    distances, indices = nbrs.kneighbors(data)

    if knn_autotune > 0:
        print('Autotuning distances')
        for j in reversed(range(N)):
            temp = sorted(distances[j])
            lMaxTempIdxs = min(knn_autotune, len(temp))
            if lMaxTempIdxs == 0 or temp[lMaxTempIdxs] == 0:
                distances[j] = 0
            else:
                distances[j] = np.divide(distances[j], temp[lMaxTempIdxs])

    # Adjacency matrix
    rows = np.zeros(N * knn, dtype=np.int32)
    cols = np.zeros(N * knn, dtype=np.int32)
    dists = np.zeros(N * knn)
    location = 0
    for i in range(N):
        inds = range(location, location + knn)
        rows[inds] = indices[i, :]
        cols[inds] = i
        dists[inds] = distances[i, :]
        location += knn
    if epsilon > 0:
        W = csr_matrix( (dists, (rows, cols)), shape=[N, N] )
    else:
        W = csr_matrix( (np.ones(dists.shape), (rows, cols)), shape=[N, N] )

    # Symmetrize W
    W = W + W.T

    if epsilon > 0:
        # Convert to affinity (with selfloops)
        rows, cols, dists = find(W)
        rows = np.append(rows, range(N))
        cols = np.append(cols, range(N))
        dists = np.append(dists/(epsilon ** 2), np.zeros(N))
        W = csr_matrix( (np.exp(-dists), (rows, cols)), shape=[N, N] )

    # Create D
    D = np.ravel(W.sum(axis = 1))
    D[D!=0] = 1/D[D!=0]

    #markov normalization
    T = csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(W)

    return T