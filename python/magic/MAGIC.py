# author: Scott Gigante <scott.gigante@yale.edu>
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
from .utils import check_int, check_positive, check_between, check_in, check_if_not
from .logging import set_logging, log_start, log_complete, log_info, log_debug

try:
    import anndata
except ImportError:
    # anndata not installed
    pass


class MAGIC(BaseEstimator):
    """MAGIC TODO [1]

    Parameters
    ----------

    # TODO

    k : int, optional, default: 15
        number of nearest neighbors on which to build kernel

    a : int, optional, default: 10
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used

    t : int, optional, default: 'auto'
        power to which the diffusion operator is powered.
        This sets the level of diffusion. If 'auto', t is selected
        according to the knee point in the Von Neumann Entropy of
        the diffusion operator

    n_pca : int, optional, default: 100
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        roughly log(n_samples) time.

    n_jobs : integer, optional, default: 1
        The number of jobs to use for the computation.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging.
        For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for
        n_jobs = -2, all CPUs but one are used

    random_state : integer or numpy.RandomState, optional, default: None
        The generator used to initialize SMACOF (metric, nonmetric) MDS
        If an integer is given, it fixes the seed
        Defaults to the global `numpy` random number generator

    verbose : `int` or `boolean`, optional (default: 1)
        If `True` or `> 0`, print status messages

    Attributes
    ----------

    X : array-like, shape=[n_samples, n_dimensions]

    X_imputed : array-like, shape=[n_samples, n_dimensions]

    graph : graphtools.BaseGraph
        The graph built on the input data

    Examples
    --------
     # TODO

    References
    ----------
    .. [1] TODO
    """

    def __init__(self, n_pca=100, t='auto', k=12, a=10, rescale=0,
                 n_jobs=1, random_state=None, verbose=1):
        self.a = a
        self.k = k
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
            if isinstance(self.graph, graphtools.LandmarkGraph):
                return self.graph.landmark_op
            else:
                return self.graph.diff_op
        else:
            raise NotFittedError("This MAGIC instance is not fitted yet. Call "
                                 "'fit' with appropriate arguments before "
                                 "using this method.")

    def _check_params(self):
        """Check MAGIC parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as t=-1, would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        check_positive(k=self.k,
                       n_pca=self.n_pca)
        check_int(k=self.k,
                  n_jobs=self.n_jobs)
        # TODO: epsilon
        check_if_not(0, check_between, v_min=1,
                     v_max=100, rescale=self.rescale)
        check_if_not(None, check_positive, check_int,
                     n_pca=self.n_pca, self.a)
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
        self.X = X

        if self.X.shape[1] <= self.n_pca:
            n_pca = None
        else:
            n_pca = self.n_pca

        if self.graph is None:

            log_start("graph and diffusion operator")
            self.graph = graphtools.Graph(
                X,
                n_pca=n_pca,
                knn=self.k + 1,
                decay=self.a,
                thresh=self.alpha_threshold,
                n_jobs=self.n_jobs,
                verbose=self.verbose,
                random_state=self.random_state)
            log_complete("graph and diffusion operator")
        else:
            # check the user hasn't changed parameters manually
            try:
                self.graph.set_params(
                    decay=self.a, knn=self.k + 1,
                    n_jobs=self.n_jobs, verbose=self.verbose, n_pca=n_pca,
                    thresh=self.alpha_threshold,
                    random_state=self.random_state)
                log_info("Using precomputed graph and diffusion operator...")
            except ValueError:
                # something changed that should have invalidated the graph
                self.graph = None
                return self.fit(X)
        return self

    def impute(self, X, t_max=100, plot_optimal_t=False, ax=None):
        # TODO: landmarks?
        log_start("imputation")
        if self.t == 'auto':
            t = self.optimal_t(X, t_max=t_max, plot=plot_optimal_t, ax=ax)
            log_info("Automatically selected t = {}".format(t))
        else:
            t = self.t
        if X.shape[1] > self.n_pca:
            pca = PCA(self.n_pca)
            X = pca.transform(X)

        for _ in range(t):
            X = self.diff_op.dot(X)

        if pca is not None:
            X = pca.inverse_transform(X)
        log_complete("imputation")
        return X

    def transform(self, X=None, t_max=100, plot_optimal_t=False, ax=None):
        """Computes the imputed values

        Parameters
        ----------
        X : array, optional, shape=[n_samples, n_features]
            input data with `n_samples` samples and `n_dimensions`
            dimensions. By default, we use the data used to call
            `MAGIC.fit()`.
            Accepted data types: `numpy.ndarray`,
            `scipy.sparse.spmatrix`, `pd.DataFrame`, `anndata.AnnData`. If
            `knn_dist` is 'precomputed', `data` should be a n_samples x
            n_samples distance or affinity matrix

        t_max : int, optional, default: 100
            maximum t to test if `t` is set to 'auto'

        plot_optimal_t : boolean, optional, default: False
            If true and `t` is set to 'auto', plot the Von Neumann
            entropy used to select t

        ax : matplotlib.axes.Axes, optional
            If given and `plot_optimal_t` is true, plot will be drawn
            on the given axis.

        Returns
        -------
        embedding : array, shape=[n_samples, n_dimensions]
        The cells embedded in a lower dimensional space using PHATE
        """
        if self.graph is None:
            raise NotFittedError("This MAGIC instance is not fitted yet. Call "
                                 "'fit' with appropriate arguments before "
                                 "using this method.")
        elif X is not None and not np.all(X == self.X):
            warnings.warn(UserWarning, "Running MAGIC.transform on different "
                          "data to that which was used for MAGIC.fit may not "
                          "produce sensible output.")
            return self.impute(
                X, t_max=t_max, plot_optimal_t=plot_optimal_t, ax=ax)
        else:
            if self.X_imputed is None:
                self.X_imputed = self.impute(
                    self.graph.data_nu, t_max=t_max,
                    plot_optimal_t=plot_optimal_t, ax=ax)
                self.X_imputed = self.graph.inverse_transform(self.X_imputed)
            else:
                log_info("Using precomputed imputation...")
            # if isinstance(self.graph, graphtools.LandmarkGraph):
            #     log_debug("Extending to original data...")
            #     return self.graph.interpolate(self.X_imputed)
            return self.X_imputed

    def fit_transform(self, X, **kwargs):
        """Computes the diffusion operator and the imputed values

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
        X_imputed = self.transform(**kwargs)
        log_complete('MAGIC')
        return X_imputed

    def optimal_t(self, data, plot=False,
                  t_max=20, compute_t_n_genes=500):
        """Find the optimal value of t

        Selects the optimal value of t based on the R squared of
        successively higher values of t

        Parameters
        ----------
        t_max : int, default: 100
            Maximum value of t to test

        plot : boolean, default: False
            If true, plots the Von Neumann Entropy and knee point

        ax : matplotlib.Axes, default: None
            If plot=True and ax is not None, plots the VNE on the given axis
            Otherwise, creates a new axis and displays the plot

        Returns
        -------
        t_opt : int
            The optimal value of t
        """
        log_start("optimal t")
        t_opt, r2_vec = compute_optimal_t(data, self.diff_op,
                                          t_max=t_max,
                                          compute_t_n_genes=compute_t_n_genes)
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
