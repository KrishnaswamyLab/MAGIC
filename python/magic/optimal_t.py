import numpy as np
from scipy import sparse


def compute_optimal_t(data, diff_op, make_plots=True, t_max=32, n_genes=None):

    if n_genes is None:
        n_genes = data.shape[1]

    if not sparse.issparse(diff_op):
        diff_op = sparse.csr_matrix(diff_op)

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

    complete = False
    for i in range(t_max):
        data_imputed = diff_op.dot(data_imputed)
        data_curr = data_imputed
        data_curr = np.divide(data_curr, np.broadcast_to(
            data_curr.sum(axis=0), data_curr.shape))
        r2 = rsquare(data_prev, data_curr)
        r2_vec[i] = 1 - r2
        if r2_vec[i] < 0.05:
            complete = True
            t_opt = i + 2
            break
        data_prev = data_curr

    if complete:
        print('optimal t = ', t_opt)
    else:
        # r2 has not converged
        t_opt = t_max
        print("Warning: optimal t > t_max ({})".format(t_opt))

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

    # rmse = np.sqrt(np.power((y - f).mean(), 2))

    return r2  # , rmse
