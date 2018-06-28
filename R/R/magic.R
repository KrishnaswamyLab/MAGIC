#' @title Perform MAGIC on a data matrix
#'
#' @param data input data matrix
#' @param genes character or integer vector, default: NULL
#' vector of column names or column indices for which to return smoothed data
#' If 'all_genes' or NULL, the entire smoothed matrix is returned
#' @param k int, optional, default: 5
#' number of nearest neighbors on which to build kernel
#' @param alpha int, optional, default: 15
#' sets decay rate of kernel tails.
#' If NULL, alpha decaying kernel is not used
#' @param t int, optional, default: 'auto'
#' power to which the diffusion operator is powered
#' sets the level of diffusion
#' @param npca number of PCA components that should be used; default: 20.
#' @param rescale_percent To which percentile should the data be re-scaled.
#' Expects an integer between 0 and 100. If 0, no rescaling is done.
#' Note: Do not set this higher than 0 if your data has negative values e.g.
#' log transformed data.
#' Default: 0.
#' @param init magic object, optional
#' object to use for initialization. Avoids recomputing
#' intermediate steps if parameters are the same.
#' @param t.max int, optional, default: 20
#' Maximum value of t to test for automatic t selection.
#' @param knn.dist.method string, optional, default: 'euclidean'.
#' recommended values: 'euclidean', 'cosine'
#' Any metric from `scipy.spatial.distance` can be used
#' distance metric for building kNN graph.
#' @param verbose `int` or `boolean`, optional (default : 1)
#' If `TRUE` or `> 0`, print verbose updates.
#' @param n.jobs `int`, optional (default: 1)
#' The number of jobs to use for the computation.
#' If -1 all CPUs are used. If 1 is given, no parallel computing code is
#' used at all, which is useful for debugging.
#' For n_jobs below -1, (n.cpus + 1 + n.jobs) are used. Thus for
#' n_jobs = -2, all CPUs but one are used
#' @param seed int or `NULL`, random state (default: `NULL`)
#' @export
magic <- function(data,
                  genes=NULL,
                  k = 5,
                  alpha = 15,
                  t = 'auto',
                  npca=100,
                  rescale_percent=0,
                  init=NULL,
                  knn.dist.method='euclidean',
                  verbose=1,
                  n.jobs=1) {
  # check installation
  if (!reticulate::py_module_available(module = "magic")) {
    install.magic()
  }
  tryCatch(pymagic, error = function(e) load_pymagic())
  k <- as.integer(k)
  t.max <- as.integer(t.max)
  n.jobs <- as.integer(n.jobs)

  if (is.numeric(n.landmark)) {
    n.landmark <- as.integer(n.landmark)
  } else if (!is.null(n.landmark) && is.na(n.landmark)) {
    n.landmark <- NULL
  }
  if (is.numeric(npca)) {
    npca <- as.integer(npca)
  } else if (!is.null(npca) && is.na(npca)) {
    npca <- NULL
  }
  if (is.numeric(alpha)) {
    alpha <- as.double(alpha)
  } else if (!is.null(alpha) && is.na(alpha)) {
    alpha <- NULL
  }
  if (is.numeric(t)) {
    t <- as.integer(t)
  } else if (is.null(t) || is.na(t)) {
    t <- 'auto'
  }
  if (is.numeric(rescale_percent)) {
    rescale_percent <- as.double(rescale_percent)
  } else if (is.null(rescale_percent) || is.na(rescale_percent)) {
    rescale_percent <- 0
  }
  if (is.numeric(seed)) {
    seed <- as.integer(seed)
  } else if (!is.null(seed) && is.na(seed)) {
    seed <- NULL
  }
  if (is.numeric(verbose)) {
    verbose <- as.integer(verbose)
  }
  if (!methods::is(data, "Matrix")) {
    data <- as.matrix(data)
  }
  if (is.numeric(genes)) {
    genes <- as.integer(genes)
    gene_names <- colnames(data)[genes]
    genes <- genes - 1
  } else if (!is.null(genes) && is.na(genes)) {
    genes <- NULL
    gene_names <- colnames(data)
  } else if (length(genes) == 1 && genes == "all_genes") {
    gene_names <- colnames(data)
  } else {
    # character vector
    gene_names <- colnames(data)[colnames(data) %in% genes]
  }

  # store parameters
  params <- list("data" = data, "k" = k, "alpha" = alpha, "t" = t,
                 "npca" = npca, "knn.dist.method" = knn.dist.method,
                 "rescale_percent" = rescale_percent)
  # use pre-initialized values if given
  operator <- NULL
  if (!is.null(init)) {
    if (!methods::is(init, "magic")) {
      warning("object passed to init is not a phate object")
    } else {
      operator <- init$operator
      operator$set_params(k = k,
                          a = alpha,
                          t = t,
                          n_pca = npca,
                          knn_dist = knn.dist.method,
                          n_jobs = n.jobs,
                          random_state = seed,
                          verbose = verbose,
                          rescale = rescale_percent)
    }
  }
  if (is.null(operator)) {
    operator <- pymagic$MAGIC(k = k,
                              a = alpha,
                              t = t,
                              n_pca = npca,
                              knn_dist = knn.dist.method,
                              n_jobs = n.jobs,
                              random_state = seed,
                              verbose = verbose,
                              rescale = rescale_percent)
  }
  result <- operator$fit_transform(data,
                                   genes = genes,
                                   t_max = t.max)
  colnames(result) <- gene_names
  rownames(result) <- rownames(data)
  result <- list("result" = result, "operator" = operator,
                 "params" = params)
  class(result) <- c("magic", "list")
  return(result)
}


#' Print a MAGIC object
#'
#' This avoids spamming the user's console with a list of many large matrices
#'
#' @param x A fitted MAGIC object
#' @param ... Arguments for print()
#' @examples
#' if (reticulate::py_module_available("magic")) {
#'
#' # data(tree.data)
#' # We use a smaller tree to make examples run faster
#' data(tree.data.small)
#' phate.tree <- phate(tree.data.small$data)
#' print(phate.tree)
#' ## MAGIC with elements
#' ## $result : (3000, 2)
#' ## $operator : Python MAGIC operator
#' ## $params : list with elements (data, k, alpha, t, n.landmark, ndim,
#' ##                               potential.method, npca, mds.method,
#' ##                               knn.dist.method, mds.dist.method)
#'
#' }
#' @rdname print
#' @method print phate
#' @export
print.phate <- function(x, ...) {
  result <- paste0("MAGIC with elements\n",
                   "  $result : (", nrow(x$embedding), ", ",
                   ncol(x$embedding), ")\n",
                   "  $operator : Python MAGIC operator\n",
                   "  $params : list with elements (",
                   paste(names(x$params), collapse = ", "), ")")
  cat(result)
}

#' Summarize a MAGIC object
#'
#' @param object A fitted MAGIC object
#' @param ... Arguments for summary()
#' @examples
#' if (reticulate::py_module_available("phate")) {
#'
#' # data(tree.data)
#' # We use a smaller tree to make examples run faster
#' data(tree.data.small)
#' phate.tree <- phate(tree.data.small$data)
#' summary(phate.tree)
#' ## MAGIC embedding
#' ## k = 5, alpha = NULL, t = 58
#' ## Data: (3000, 100)
#' ## Embedding: (3000, 2)
#'
#' }
#' @rdname summary
#' @method summary phate
#' @export
summary.phate <- function(object, ...) {
  summary(object$result)
}

#' Convert a MAGIC object to a matrix
#'
#' Returns the smoothed data matrix
#'
#' @param x A fitted MAGIC object
#' @param ... Arguments for as.matrix()
#' @rdname as.matrix
#' @method as.matrix phate
#' @export
as.matrix.phate <- function(x, ...) {
  x$result
}
#' Convert a MAGIC object to a data.frame
#'
#' Returns the smoothed data matrix
#'
#' @param x A fitted MAGIC object
#' @param ... Arguments for as.data.frame()
#' @rdname as.data.frame
#' @method as.data.frame phate
#' @export
as.data.frame.phate <- function(x, ...) {
  as.data.frame(as.matrix(x), ...)
}

