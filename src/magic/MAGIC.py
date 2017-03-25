import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix, find
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.manifold.t_sne import _joint_probabilities, _joint_probabilities_nn
from sklearn.decomposition import PCA
from scipy.spatial.distance import squareform
from sklearn.neighbors import NearestNeighbors

def magic(data, n_pca_components=20, random_pca=True, 
          t=6, knn=30, knn_autotune=10, epsilon=1, rescale=99):

    if n_pca_components != None:
        print('doing PCA')
        pca_projected_data = run_pca(data, n_components=n_pca_components, random=random_pca)
    else:
        pca_projected_data = data

    #run diffusion maps to get markov matrix
    L = compute_markov(pca_projected_data, knn=knn, epsilon=epsilon, 
                       distance_metric='euclidean', knn_autotune=knn_autotune)

    #remove tsne kernel for now
    # else:
    #     distances = pairwise_distances(pca_projected_data, squared=True)
    #     if k_knn > 0:
    #         neighbors_nn = np.argsort(distances, axis=0)[:, :k_knn]
    #         P = _joint_probabilities_nn(distances, neighbors_nn, perplexity, 1)
    #     else:
    #         P = _joint_probabilities(distances, perplexity, 1)
    #     P = squareform(P)

    #     #markov normalize P
    #     L = np.divide(P, np.sum(P, axis=1))

    #get imputed data matrix -- by default use data_norm but give user option to pick
    new_data, L_t = impute_fast(data, L, t, rescale_percent=rescale)

    return new_data


def run_pca(data, n_components=100, random=True):

    solver = 'randomized'
    if random != True:
        solver = 'full'

    pca = PCA(n_components=n_components, svd_solver=solver)
    return pca.fit_transform(data)


def impute_fast(data, L, t, rescale_percent=0, L_t=None, tprev=None):

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

    #rescale data
    if rescale_percent != 0:
        if len(np.where(data_new < 0)[0]) > 0:
            print('Rescaling should not be performed on log-transformed '
                  '(or other negative) values. Imputed data returned unscaled.')
            return data_new, L_t
            
        M99 = np.percentile(data, rescale_percent, axis=0)
        M100 = data.max(axis=0)
        indices = np.where(M99 == 0)[0]
        M99[indices] = M100[indices]
        M99_new = np.percentile(data_new, rescale_percent, axis=0)
        M100_new = data_new.max(axis=0)
        indices = np.where(M99_new == 0)[0]
        M99_new[indices] = M100_new[indices]
        max_ratio = np.divide(M99, M99_new)
        data_new = np.multiply(data_new, np.tile(max_ratio, (len(data), 1)))
    
    return data_new, L_t


def compute_markov(data, knn=10, epsilon=1, distance_metric='euclidean', knn_autotune=0):

    N = data.shape[0]

    # Nearest neighbors
    print('Computing distances')
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
    print('Computing kernel')
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


def optimal_t(data, th=0.001):
    S = np.linalg.svd(data)
    S = np.power(S, 2)
    nse = np.zeros(32)

    for t in range(32):
        S_t = np.power(S, t)
        P = np.divide(S_t, np.sum(S_t, axis=0))
        nse[t] = np.sum(P[np.where(P > th)[0]])

