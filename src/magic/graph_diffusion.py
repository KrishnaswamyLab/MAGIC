import numpy as np
import pandas as pd
import time
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix, find, issparse
from scipy.sparse.linalg import eigs
from numpy.linalg import norm
from magic.bimarkov import bimarkov
from magic.GetEigs import GetEigs

def run_diffusion_map(data, knn=10, normalization='smarkov', 
                      epsilon=1, n_diffusion_components=10, 
                      distance_metric='minkowski', knn_autotune=0):
    """ Run diffusion maps on the data. This implementation is based on the 
        diffusion geometry library in Matlab: https://services.math.duke.edu/~mauro/code.html#DiffusionGeom
    :param data: Data matrix of samples X features
    :param knn: Number of neighbors for graph construction to determine distances between cells
    :param normalization: method for normalizing the matrix of weights
         'bimarkov'            force row and column sums to be 1
         'markov'              force row sums to be 1
         'smarkov'             symmetric conjugate to markov
         'beltrami'            Laplace-Beltrami normalization ala Coifman-Lafon
         'sbeltrami'           symmetric conjugate to beltrami
         'FokkerPlanck'        Fokker-Planck normalization
         'sFokkerPlanck'       symmetric conjugate to Fokker-Planck normalization
    :param epsilon: Gaussian standard deviation for converting distances to affinities
    :param n_diffusion_components: Number of diffusion components to generate
    :return: Dictionary containing diffusion operator, weight matrix, 
             diffusion eigen vectors, and diffusion eigen values
    """

    if normalization not in ['bimarkov', 'smarkov', 'markov', 'sbeltrami', 'beltrami',
        'FokkerPlanck', 'sFokkerPlanck']:
        raise RuntimeError('Unsupported normalization. Please refer to the docstring for the supported methods')

    # Log
    print('Running Diffusion maps with the following parameters:')
    print('Normalization: %s' % normalization)

    start = time.process_time()
    N = data.shape[0]
    #Check if sparse, square matrix was input and treat as W
    if issparse(data) and data.shape[0] == data.shape[1]:
        print('Using precomputed affinity matrix')
        W = data
        
    else:
        print('Number of nearest neighbors k: %d' % knn)
        print('Epsilon: %.4f' % epsilon)

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

    P = None

    #Go through the various normalizations
    if normalization == 'bimarkov':
        print('(bimarkov) ... ')
        T = bimarkov(W)

    elif normalization == 'smarkov':
        print('(symmetric markov) ... ')

        D = csr_matrix((np.sqrt(D), (range(N), range(N))),  shape=[N, N])
        P = D
        T = D.dot(W).dot(D)

        T = (T + T.T) / 2
    
    elif normalization == 'markov':
        print('(markov) ... ')

        T = csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(W)
    
    elif normalization == 'sbeltrami':
        print('(symmetric beltrami) ... ')
    
        P = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        K = P.dot(W).dot(P)
    
        D = np.ravel(K.sum(axis = 1))
        D[D!=0] = 1/D[D!=0]
    
        D = csr_matrix((D, (range(N), range(N))),  shape=[N, N])
        P = D
        T = D.dot(K).dot(D)

        T = (T + T.T) / 2    # iron out numerical wrinkles
    
    elif normalization == 'beltrami':
        print('(beltrami) ... ')

        D = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        K = D.dot(W).dot(D)
    
        D = np.ravel(K.sum(axis = 1))
        D[D!=0] = 1/D[D!=0]
    
        V = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        T = V.dot(K)
    
    elif normalization == 'FokkerPlanck':
        print('(FokkerPlanck) ... ')
    
        D = csr_matrix((np.sqrt(D), (range(N), range(N))),  shape=[N, N])
        K = D.dot(W).dot(D)
    
        D = np.ravel(K.sum(axis = 1))
        D[D!=0] = 1/D[D!=0]
    
        D = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        T = D.dot(K)
    
    elif normalization == 'sFokkerPlanck':
        print('(sFokkerPlanck) ... ')

        D = csr_matrix((np.sqrt(D), (range(N), range(N))),  shape=[N, N])
        K = D.dot(W).dot(D)
    
        D = np.ravel(K.sum(axis = 1))
        D[D!=0] = 1/D[D!=0]
    
        D = csr_matrix((np.sqrt(D), (range(N), range(N))),  shape=[N, N])
        P = D
        T = D.dot(K).dot(D)
    
        T = (T + T.T) / 2

    else:
        raise RuntimeError("unknown normalization")

    if normalization != 'bimarkov':
        print('%.2f seconds' % (time.process_time()-start))

    # Eigen value decomposition
    V, D = GetEigs(T, n_diffusion_components, P, take_diagonal=1)

    return {'T': T, 'W': W, 'EigenVectors': V, 'EigenValues': D}