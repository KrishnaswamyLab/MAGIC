"""
Markov Affinity-based Graph Imputation of Cells (MAGIC)
"""

# author: 
# (C) 2017 Krishnaswamy Lab GPLv2

from __future__ import print_function, division, absolute_import

import numpy as np
import graphtools
from sklearn.base import BaseEstimator
from sklearn.exceptions import NotFittedError
import warnings
import preprocessing.library_size_normalize as normalize
import matplotlib.pyplot as plt
from scipy import sparse, stats


from .utils import check_int, check_positive, check_between, check_in, check_if_not, convert_to_same_format
from .logging import set_logging, log_start, log_complete, log_info, log_debug

try:
    import anndata
except ImportError:
    # anndata not installed
    pass


class MAGIC(BaseEstimator):
    """MAGIC operator which performs dimensionality reduction.

    Markov Affinity-based Graph Imputation of Cells (MAGIC) is an 
    algorithm for denoising and transcript recover of single cells
    applied to single-cell RNA sequencing data.

    Parameters
    ----------

    k : int, optional, default: 10
        number of nearest neighbors on which to build kernel

    a : int, optional, default: 10
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used

    t : int, optional, default: 'auto'
        power to which the diffusion operator is powered.
        This sets the level of diffusion. If 'auto', t is selected
        according to the R squared of the diffused data

    n_pca : int, optional, default: 100
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        roughly log(n_samples) time.

    knn_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean', 'cosine', 'precomputed'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for building kNN graph. If 'precomputed',
        `data` should be an n_samples x n_samples distance or
        affinity matrix

    n_jobs : integer, optional, default: 1
        The number of jobs to use for the computation.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging.
        For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for
        n_jobs = -2, all CPUs but one are used

    random_state : integer or numpy.RandomState, optional, default: None
        The generator used to initialize random PCA
        If an integer is given, it fixes the seed
        Defaults to the global `numpy` random number generator

    magic_type : string, optional, default: 'strict'
        The type of magic implementation used
        'strict' ensures that the algorithm specified
        in the Magic paper is used. (W^t * D)
        'fast' is recommended when the dimensions of 
        the Markov transition matrix W (to be powered
        by t) are smaller than those of the preprocessed
        data matrix D. (W^t) * D

    verbose : `int` or `boolean`, optional (default: 1)
        If `True` or `> 0`, print status messages

    Attributes
    ----------

    X : array-like, shape=[n_samples, n_dimensions]
        Input data

    X_magic : array-like, shape=[n_samples, n_dimensions]
        Output data

    graph : graphtools.BaseGraph
        The graph built on the input data

    Examples
    --------
    >>> # TODO: better example data
    >>> import magic
    >>> tree_data, tree_clusters = phate.tree.gen_dla(n_dim=100,
                                                      n_branch=20,
                                                      branch_length=100)
    >>> tree_data.shape
    (2000, 100)
    >>> magic_operator = magic.MAGIC(k=5, a=20, t=150)
    >>> tree_magic = phate_operator.fit_transform(tree_data)
    >>> tree_magic.shape
    (2000, 100)
    >>> import phate
    >>> import matplotlib.pyplot as plt
    >>> phate_operator = phate.PHATE(knn_dist='precomputed')
    >>> tree_phate = phate_operator.fit_transform(magic_operator.graph.kernel)
    >>> # plt.scatter(tree_phate[:,0], tree_phate[:,1], c=tree_magic[:,0])
    >>> # plt.show()

  
    """

    def __init__(self, k=10, a=10, t='auto', n_pca=100, knn_dist='euclidean',
                 n_jobs=1, random_state=None, magic_type='strict', verbose=1):
        self.k = k
        self.a = a
        self.t = t
        self.n_pca = n_pca
        self.knn_dist = knn_dist
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.magic_type = magic_type

        self.graph = None
        self.X = None
        self._check_params()
        self.verbose = verbose
        set_logging(verbose)


    @property
    def diff_op(self):
        """The diffusion operator calculated from the data
        """
        if self.graph is not None:
            return self.graph.diff_op
        else:
            raise NotFittedError("This MAGIC instance is not fitted yet. Call "
                                 "'fit' with appropriate arguments before "
                                 "using this method.")


    def _check_params(self):
        """Check MAGIC parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as mds='mmds', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        check_positive(k=self.k,
                       a=self.a,
                       n_pca=self.n_pca)
        check_int(k=self.k,
                  n_jobs=self.n_jobs)
        # TODO: epsilon
        check_between(v_min=0,
                      v_max=100,
                      rescale=self.rescale)
        check_if_not(None, check_positive, check_int,
                     n_pca=self.n_pca)
        check_if_not('auto', check_positive, check_int,
                     t=self.t)
        check_in(['euclidean', 'precomputed', 'cosine', 'correlation',
                  'cityblock', 'l1', 'l2', 'manhattan', 'braycurtis',
                  'canberra', 'chebyshev', 'dice', 'hamming', 'jaccard',
                  'kulsinski', 'mahalanobis', 'matching', 'minkowski',
                  'rogerstanimoto', 'russellrao', 'seuclidean',
                  'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'],
                 knn_dist=self.knn_dist)


    def fit(self, X):
        """Computes the diffusion operator

        Parameters
        ----------
        X : array, shape=[n_samples, n_features]
            input data with `n_samples` samples and `n_dimensions`
            dimensions. Accepted data types: `numpy.ndarray`,
            `scipy.sparse.spmatrix`, `pd.DataFrame`, `anndata.AnnData`. If
            `knn_dist` is 'precomputed', `data` should be a n_samples x
            n_samples distance or affinity matrix

        Returns
        -------
        magic_operator : MAGIC
            The estimator object
        """
        try:
            if isinstance(X, anndata.AnnData):
                X = X.X
        except NameError:
            # anndata not installed
            pass

        if self.knn_dist == 'precomputed':
            if isinstance(X, sparse.coo_matrix):
                X = X.tocsr()
            if X[0, 0] == 0:
                precomputed = "distance"
            else:
                precomputed = "affinity"
            log_info("Using precomputed {} matrix...".format(precomputed))
            n_pca = None
        else:
            precomputed = None
            if X.shape[1] <= self.n_pca:
                n_pca = None
            else:
                n_pca = self.n_pca

        if self.rescale != 0 and np.any(X < 0):
            warnings.warn(
                'Rescaling should not be performed on log-transformed (or'
                ' other negative) values. Setting `rescale=0`.',
                RuntimeWarning)
            self.rescale = 0
        if self.graph is not None:
            if self.X is not None and not (X != self.X).sum() == 0:
                """
                If the same data is used, we can reuse existing kernel and
                diffusion matrices. Otherwise we have to recompute.
                """
                self.graph = None
            else:
                try:
                    self.graph.set_params(
                        decay=self.a, knn=self.k + 1, distance=self.knn_dist,
                        precomputed=precomputed,
                        n_jobs=self.n_jobs, verbose=self.verbose, n_pca=n_pca,
                        thresh=1e-4, random_state=self.random_state)
                    log_info(
                        "Using precomputed graph and diffusion operator...")
                except ValueError as e:
                    # something changed that should have invalidated the graph
                    log_debug("Reset graph due to {}".format(str(e)))
                    self.graph = None

        self.X = X

        if self.graph is None:
            log_start("graph and diffusion operator")
            self.graph = graphtools.Graph(
                X,
                n_pca=n_pca,
                knn=self.k + 1,
                decay=self.a,
                thresh=1e-4,
                n_jobs=self.n_jobs,
                verbose=self.verbose,
                random_state=self.random_state)
            log_complete("graph and diffusion operator")

        return self


    def transform(self, X=None, t_max=20, plot_optimal_t=False, ax=None):
        """Computes the position of the cells in the embedding space

        Parameters
        ----------
        X : array, optional, shape=[n_samples, n_features]
            input data with `n_samples` samples and `n_dimensions`
            dimensions. Not required, since PHATE does not currently embed
            cells not given in the input matrix to `PHATE.fit()`.
            Accepted data types: `numpy.ndarray`,
            `scipy.sparse.spmatrix`, `pd.DataFrame`, `anndata.AnnData`. If
            `knn_dist` is 'precomputed', `data` should be a n_samples x
            n_samples distance or affinity matrix

        t_max : int, optional, default: 100
            maximum t to test if `t` is set to 'auto'

        plot_optimal_t : boolean, optional, default: False
            If true and `t` is set to 'auto', plot the R squared used to
            select t

        ax : matplotlib.axes.Axes, optional
            If given and `plot_optimal_t` is true, plot will be drawn
            on the given axis.

        Returns
        -------
        embedding : array, shape=[n_samples, n_dimensions]
        The cells embedded in a lower dimensional space using PHATE
        """
        if self.graph is None:
            raise NotFittedError("This PHATE instance is not fitted yet. Call "
                                 "'fit' with appropriate arguments before "
                                 "using this method.")
        elif X is not None and np.sum(X != self.X) > 0:
            warnings.warn(UserWarning, "Running MAGIC.transform on different "
                          "data to that which was used for MAGIC.fit may not "
                          "produce sensible output, unless it comes from the "
                          "same manifold.")
            X_magic = self.impute(X)
            X_magic = convert_to_same_format(X_magic, X)
            return X_magic
        else:
            self.X_magic = self.impute(self.graph, t_max=t_max,
                                       plot=plot_optimal_t, ax=ax)
            self.X_magic = self.rescale_data(self.X_magic, self.graph.data)
            self.X_magic = convert_to_same_format(self.X_magic, self.X)
            return self.X_magic


    def fit_transform(self, X, **kwargs):
        """Computes the diffusion operator and the position of the cells in the
        embedding space

        Parameters
        ----------
        X : array, shape=[n_samples, n_features]
            input data with `n_samples` samples and `n_dimensions`
            dimensions. Accepted data types: `numpy.ndarray`,
            `scipy.sparse.spmatrix`, `pd.DataFrame`, `anndata.AnnData` If
            `knn_dist` is 'precomputed', `data` should be a n_samples x
            n_samples distance or affinity matrix

        kwargs : further arguments for `PHATE.transform()`
            Keyword arguments as specified in :func:`~phate.PHATE.transform`

        Returns
        -------
        embedding : array, shape=[n_samples, n_dimensions]
            The cells embedded in a lower dimensional space using PHATE
        """
        log_start('MAGIC')
        self.fit(X)
        X_magic = self.transform(**kwargs)
        log_complete('MAGIC')
        return X_magic


    def rsquare(self, data, data_prev):
        """
        Returns
        -------

        r2 : R squared value

        data_curr : transformed data for next time
        """
        data_curr = np.divide(data, np.broadcast_to(
            data.sum(axis=0), data.shape)).reshape(-1)
        data_curr = data_curr - data_curr.min()
        if data_prev is not None:
            r = stats.linregress(data_prev, data_curr).rvalue
            r2 = 1 - r**2
        else:
            r2 = None
        return r2, data_curr


    def impute(self, data, t_max=20, plot=False, ax=None):
        """Impute with PCA

        Parameters
        ----------
        data : graphtools.Graph, graphtools.Data or array-like
        """
        if not isinstance(data, graphtools.base.Data):
            data = graphtools.base.Data(data, n_pca=self.n_pca)
        data_imputed = data.data_nu

        if self.t == 'auto':
            _, data_prev = self.rsquare(data_imputed, data_prev=None)
            r2_vec = []
            t_opt = None
        else:
            t_opt = self.t

        log_start("imputation")
        i = 0
        while (t_opt is None and i < t_max) or (i < t_opt):
            i += 1
            data_imputed = self.diff_op.dot(data_imputed)
            if self.t == 'auto':
                r2, data_prev = self.rsquare(data_imputed, data_prev)
                r2_vec.append(r2)
                log_debug("{}: {}".format(i, r2_vec))
                if r2 < 0.05 and t_opt is None:
                    t_opt = i + 2
                    log_info("Automatically selected t = {}".format(t_opt))
        log_complete("imputation")

        if plot:
            # continue to t_max
            log_start("optimal t plot")
            if t_opt is None:
                # never converged
                warnings.warn("optimal t > t_max ({})".format(t_max),
                              RuntimeWarning)
            else:
                data_overimputed = data_imputed
                while i < t_max:
                    i += 1
                    data_overimputed = self.diff_op.dot(data_overimputed)
                    r2, data_prev = self.rsquare(data_overimputed, data_prev)
                    r2_vec.append(r2)

            # create axis
            if ax is None:
                fig, ax = plt.subplots()
                show = True
            else:
                show = False

            # plot
            x = np.arange(len(r2_vec)) + 1
            ax.plot(x, r2_vec)
            ax.plot(t_opt, r2_vec[t_opt - 1], 'ro', markersize=10,)
            ax.plot(x, np.full(len(r2_vec), 0.05), 'k--')
            ax.set_xlabel('t')
            ax.set_ylabel('1 - R^2(data_{t},data_{t-1})')
            ax.set_xlim([1, len(r2_vec)])
            ax.set_ylim([0, 1])
            plt.tight_layout()
            log_complete("optimal t plot")
            if show:
                plt.show()

        data_imputed = data.inverse_transform(data_imputed)
        return data_imputed

    