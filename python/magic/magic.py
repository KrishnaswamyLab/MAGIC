"""
Potential of Heat-diffusion for Affinity-based Trajectory Embedding (PHATE)
"""

# author: Daniel Burkhardt <daniel.burkhardt@yale.edu>
# (C) 2017 Krishnaswamy Lab GPLv2
from __future__ import print_function, division, absolute_import

import numpy as np
import graphtools
from sklearn.base import BaseEstimator
from sklearn.exceptions import NotFittedError
from sklearn.decomposition import PCA
import warnings
import matplotlib.pyplot as plt

from .optimal_t import compute_optimal_t
from .utils import check_int, check_positive, check_between, check_in, check_if_not, convert_to_same_format
from .logging import set_logging, log_start, log_complete, log_info, log_debug

try:
    import anndata
except ImportError:
    # anndata not installed
    pass


class MAGIC(BaseEstimator):
    """MAGIC operator which performs dimensionality reduction.

    TODO: description and citation

    Potential of Heat-diffusion for Affinity-based Trajectory Embedding
    (PHATE) embeds high dimensional single-cell data into two or three
    dimensions for visualization of biological progressions as described
    in Moon et al, 2017 [1]_.

    Parameters
    ----------

    k : int, optional, default: 10
        number of nearest neighbors on which to build kernel

    a : int, optional, default: 10
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used

    rescale : integer between 0 and 100, optional, default: 99
        Percentile to rescale data to after running MAGIC
        such that the output data has the same range and the
        input data. If 0, no rescaling is performed

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

    References
    ----------
    .. [1] TODO: Update reference Moon KR, van Dijk D, Zheng W, *et al.* (2017),
        *PHATE: A Dimensionality Reduction Method for Visualizing Trajectory
        Structures in High-Dimensional Biological Data*,
        `BioRxiv <http://biorxiv.org/content/early/2017/03/24/120378>`_.
    """

    def __init__(self, k=10, a=10, rescale=99,
                 t='auto', n_pca=100,
                 n_jobs=1, random_state=None, verbose=1):
        self.k = k
        self.a = a
        self.t = t
        self.n_pca = n_pca
        self.n_jobs = n_jobs
        self.random_state = random_state

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
        """Check PHATE parameters

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
        check_if_not(0, check_between, v_min=1,
                     v_max=100, rescale=self.rescale)
        check_if_not(None, check_positive, check_int,
                     n_pca=self.n_pca)
        check_if_not('auto', check_positive, check_int,
                     t=self.t)

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
        if self.X is not None and not np.all(X == self.X):
            """
            If the same data is used, we can reuse existing kernel and
            diffusion matrices. Otherwise we have to recompute.
            """
            self.graph = None

        if self.knn_dist == 'precomputed':
            if isinstance(X, sparse.coo_matrix):
                X = X.tocsr()
            if X[0, 0] == 0:
                precomputed = "distance"
            else:
                precomputed = "affinity"
            n_pca = None
        else:
            precomputed = None
            if X.shape[1] <= self.n_pca:
                n_pca = None
            else:
                n_pca = self.n_pca
        if self.n_landmark is None or X.shape[0] <= self.n_landmark:
            n_landmark = None
        else:
            n_landmark = self.n_landmark

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
                    log_info("Using precomputed graph and diffusion operator...")
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

    def transform(self, X=None, t_max=100, plot_optimal_t=False, ax=None):
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
            If true and `t` is set to 'auto', plot the R squared used to select t

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
            self.X_magic = self.impute(self.graph)
            self.X_magic = self.rescale(self.X_magic, self.graph.data)
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

    def impute(self, X, t_max=20, plot=False):
        """Impute with PCA

        Parameters
        ----------
        X : graphtools.Graph or array-like
        """
        if not isinstance(X, graphtools.base.Data):
            X = graphtools.base.Data(X, n_pca=self.n_pca)
        data = X.data_nu
        if self.t == 'auto':
            t = self.optimal_t(data, t_max=t_max, plot=plot_optimal_t, ax=ax)
            log_info("Automatically selected t = {}".format(t))
        else:
            t = self.t
        for _ in range(t):
            data = self.diff_op.dot(data)
        data = X.inverse_transform(data)
        return data

    def rescale(self, data, target_data):
        if self.rescale == 0:
            return data
        elif np.any(data < 0):
            warnings.warn('Rescaling should not be performed on log-transformed (or'
                          ' other negative) values. Imputed data returned unscaled.',
                          RuntimeWarning)
            return data
        else:
            M99 = np.percentile(target_data, self.rescale, axis=0)
            M100 = target_data.max(axis=0)
            indices = np.where(M99 == 0)[0]
            M99[indices] = M100[indices]
            M99_new = np.percentile(data, self.rescale, axis=0)
            M100_new = data.max(axis=0)
            indices = np.where(M99_new == 0)[0]
            M99_new[indices] = M100_new[indices]
            max_ratio = np.divide(M99, M99_new)
            data = np.multiply(data, np.tile(max_ratio,
                                             (target_data.shape[0], 1)))
        return data

    def optimal_t(self, data, plot=False,
                  t_max=20, compute_t_n_genes=500):
        """Find the optimal value of t

        Selects the optimal value of t based on the R squared of the transformed data

        Parameters
        ----------
        t_max : int, default: 20
            Maximum value of t to test

        plot : boolean, default: False
            If true, plots the R squared and convergence point

        ax : matplotlib.Axes, default: None
            If plot=True and ax is not None, plots the R squared on the given axis
            Otherwise, creates a new axis and displays the plot

        Returns
        -------
        t_opt : int
            The optimal value of t
        """
        log_start("optimal t")
        t_opt, r2_vec = compute_optimal_t(data, self.diff_op, compute_t_make_plots=True,
                                          t_max=20, compute_t_n_genes=500)
        log_complete("optimal t")

        if plot:
            plt.figure()
            plt.plot(np.arange(1, t_max + 1), r2_vec)
            plt.plot(t_opt, r2_vec[t_opt - 1], 'ro', markersize=10,)
            plt.plot(np.arange(1, t_max + 1),
                     np.full(len(r2_vec), 0.05), 'k--')
            plt.xlabel('t')
            plt.ylabel('1 - R^2(data_{t},data_{t-1})')
            plt.xlim([1, t_max])
            plt.ylim([0, 1])
            plt.tight_layout()
    #     legend({'y' 'optimal t' 'y=0.05'});

        return t_opt
