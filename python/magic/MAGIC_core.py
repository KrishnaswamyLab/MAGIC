import numpy as np
from scipy.sparse import issparse, csr_matrix, find
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt


def magic(data, diffusion_operator=None, n_pca_components=100, random_pca=True,
          t=None, k=12, ka=4, epsilon=1, rescale=0, compute_t_make_plots=True,
          t_max=12, compute_t_n_genes=500):

    if diffusion_operator is None:
        if n_pca_components is not None:
            print('doing PCA')
            pca_projected_data = run_pca(
                data, n_components=n_pca_components, random=random_pca)
        else:
            pca_projected_data = data

        # run diffusion maps to get markov matrix
        W = compute_markov(pca_projected_data, k=k, epsilon=epsilon,
                           distance_metric='euclidean', ka=ka)
    else:
        W = diffusion_operator

    if t is None:
        t, _ = compute_optimal_t(
            data, W, make_plots=compute_t_make_plots, t_max=t_max,
            n_genes=compute_t_n_genes)

    # get imputed data matrix
    new_data, _ = impute_fast(data, W, t, rescale_percent=rescale)

    return new_data


def run_pca(data, n_components=100, random=True):

    solver = 'randomized'
    if not random:
        solver = 'full'

    pca = PCA(n_components=n_components, svd_solver=solver)
    return pca.fit_transform(data)


def impute_fast(data, W, t, rescale_percent=0, W_t=None, tprev=None):

    # convert L to full matrix
    if issparse(W):
        W = W.todense()

    # L^t
    print('MAGIC: W_t = W^t')
    if W_t is None:
        W_t = np.linalg.matrix_power(W, t)
    else:
        W_t = np.dot(W_t, np.linalg.matrix_power(W, t - tprev))

    print('MAGIC: data_new = W_t * data')
    if t > 0:
        data_new = np.array(np.dot(W_t, data))

    # rescale data
    if rescale_percent != 0:
        if len(np.where(data_new < 0)[0]) > 0:
            print('Rescaling should not be performed on log-transformed (or'
                  ' other negative) values. Imputed data returned unscaled.')
            return data_new, W_t

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

    return data_new, W_t


def compute_markov(data, k=10, epsilon=1, distance_metric='euclidean', ka=0):

    N = data.shape[0]

    # Nearest neighbors
    print('Computing distances')
    nbrs = NearestNeighbors(n_neighbors=k, metric=distance_metric).fit(data)
    distances, indices = nbrs.kneighbors(data)

    if ka > 0:
        print('Autotuning distances')
        for j in reversed(range(N)):
            temp = sorted(distances[j])
            lMaxTempIdxs = min(ka, len(temp))
            if lMaxTempIdxs == 0 or temp[lMaxTempIdxs] == 0:
                distances[j] = 0
            else:
                distances[j] = np.divide(distances[j], temp[lMaxTempIdxs])

    # Adjacency matrix
    print('Computing kernel')
    rows = np.zeros(N * k, dtype=np.int32)
    cols = np.zeros(N * k, dtype=np.int32)
    dists = np.zeros(N * k)
    location = 0
    for i in range(N):
        inds = range(location, location + k)
        rows[inds] = indices[i, :]
        cols[inds] = i
        dists[inds] = distances[i, :]
        location += k
    if epsilon > 0:
        W = csr_matrix((dists, (rows, cols)), shape=[N, N])
    else:
        W = csr_matrix((np.ones(dists.shape), (rows, cols)), shape=[N, N])

    # Symmetrize W
    W = W + W.T

    if epsilon > 0:
        # Convert to affinity (with selfloops)
        rows, cols, dists = find(W)
        rows = np.append(rows, range(N))
        cols = np.append(cols, range(N))
        dists = np.append(dists / (epsilon ** 2), np.zeros(N))
        W = csr_matrix((np.exp(-dists), (rows, cols)), shape=[N, N])

    # Create D
    D = np.ravel(W.sum(axis=1))
    D[D != 0] = 1 / D[D != 0]

    # markov normalization
    T = csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(W)

    return T


def compute_optimal_t(data, diff_op, make_plots=True, t_max=32, n_genes=None):

    if n_genes is None:
        n_genes = data.shape[1]

    if not issparse(diff_op):
        diff_op = csr_matrix(diff_op)

    if n_genes < data.shape[1]:
        idx_genes = np.random.randint(data.shape[1], size=n_genes)
    else:
        idx_genes = np.arange(data.shape[1])
    data_imputed = data
    data_imputed = data_imputed[:, idx_genes]

    if data_imputed.min() < 0:
        print('data has negative values, shifting to positive')
        data_imputed = data_imputed - data_imputed.min()

    r2_vec = np.full(t_max, np.nan)
    data_prev = data_imputed
    data_prev = np.divide(data_prev, np.broadcast_to(
        data_prev.sum(axis=0), data_prev.shape))

    print('computing optimal t')
    for i in range(t_max):
        data_imputed = diff_op.dot(data_imputed)
        data_curr = data_imputed
        data_curr = np.divide(data_curr, np.broadcast_to(
            data_curr.sum(axis=0), data_curr.shape))
        r2, rmse = rsquare(data_prev, data_curr)
        r2_vec[i] = 1 - r2
        data_prev = data_curr

    try:
        t_opt = np.min(np.where(r2_vec < 0.05)) + 2
        print('optimal t = ', t_opt)
    except ValueError:
        # r2 has not converged
        t_opt = t_max
        print("Warning: optimal t > t_max ({})".format(t_opt))

    # plot
    if make_plots:
        plt.figure()
        plt.plot(np.arange(1, t_max + 1), r2_vec)
        plt.plot(t_opt, r2_vec[t_opt - 1], 'ro', markersize=10,)
        plt.plot(np.arange(1, t_max + 1), np.full(len(r2_vec), 0.05), 'k--')
        plt.xlabel('t')
        plt.ylabel('1 - R^2(data_{t},data_{t-1})')
        plt.xlim([1, t_max])
        plt.ylim([0, 1])
        plt.tight_layout()
    #     legend({'y' 'optimal t' 'y=0.05'});

    return t_opt, r2_vec


def rsquare(y, f, c=True):

    if y.shape != f.shape:
        raise RuntimeError('Y and F must be the same size.')

    # flatten arrays
    y = np.ravel(y)
    f = np.ravel(f)

    # remove NaNs
    temp = np.intersect1d(np.where(np.isfinite(y))[
                          0], np.where(np.isfinite(f))[0])
    y = y[temp]
    f = f[temp]

    if c:
        r2 = max(0, 1 - np.power(y - f, 2).sum() /
                 np.power(y - y.mean(), 2).sum())
    else:
        r2 = 1 - np.power(y - f, 2).sum() / np.power(y, 2).sum()
        if r2 < 0:
            # http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
            print('Warning: Consider adding a constant term to your model')
            r2 = 0

    rmse = np.sqrt(np.power((y - f).mean(), 2))

    return r2, rmse
