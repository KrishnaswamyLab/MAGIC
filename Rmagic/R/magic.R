#' Perform MAGIC on a data matrix
#'
#' Markov Affinity-based Graph Imputation of Cells (MAGIC) is an
#' algorithm for denoising and transcript recover of single cells
#' applied to single-cell RNA sequencing data, as described in
#' van Dijk et al, 2018.
#'
#' @param data input data matrix
#' @param genes character or integer vector, default: NULL
#' vector of column names or column indices for which to return smoothed data
#' If 'all_genes' or NULL, the entire smoothed matrix is returned
#' @param k int, optional, default: 10
#' number of nearest neighbors on which to build kernel
#' @param alpha int, optional, default: 15
#' sets decay rate of kernel tails.
#' If NULL, alpha decaying kernel is not used
#' @param t int, optional, default: 'auto'
#' power to which the diffusion operator is powered
#' sets the level of diffusion. If 'auto', t is selected according to the
#' Procrustes disparity of the diffused data.'
#' @param npca number of PCA components that should be used; default: 100.
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
#'
#' @examples
#' if (reticulate::py_module_available("magic")) {
#'
#' data(magic_testdata)
#'
#' # Run MAGIC
#' data_magic <- magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"))
#' summary(data_magic)
#' ##       CDH1             VIM             ZEB1
#' ## Min.   :0.4303   Min.   :3.854   Min.   :0.01111
#' ## 1st Qu.:0.4444   1st Qu.:3.947   1st Qu.:0.01145
#' ## Median :0.4462   Median :3.964   Median :0.01153
#' ## Mean   :0.4461   Mean   :3.965   Mean   :0.01152
#' ## 3rd Qu.:0.4478   3rd Qu.:3.982   3rd Qu.:0.01160
#' ## Max.   :0.4585   Max.   :4.127   Max.   :0.01201
#'
#' # Plot the result with ggplot2
#' if (require(ggplot2)) {
#'   ggplot(data_magic) +
#'     geom_point(aes(x=VIM, y=CDH1, color=ZEB1))
#' }
#'
#' # Run MAGIC again returning all genes
#' # We use the last run as initialization
#' data_magic <- magic(magic_testdata, genes="all_genes", init=data_magic)
#' # Extract the smoothed data matrix to use in downstream analysis
#' data_smooth <- as.matrix(data_magic)
#'
#' }
#'
#' @export
magic <- function(data,
                  genes=NULL,
                  k = 10,
                  alpha = 15,
                  t = 'auto',
                  npca=100,
                  init=NULL,
                  t.max=20,
                  knn.dist.method='euclidean',
                  verbose=1,
                  n.jobs=1,
                  seed=NULL) {
  # check installation
  if (!reticulate::py_module_available(module = "magic")) {
    load_pymagic()
  }
  tryCatch(pymagic, error = function(e) load_pymagic())
  k <- as.integer(k)
  t.max <- as.integer(t.max)
  n.jobs <- as.integer(n.jobs)

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
  if (is.null(genes) || is.na(genes)) {
    genes <- NULL
    gene_names <- colnames(data)
  } else if (is.numeric(genes)) {
    gene_names <- colnames(data)[genes]
    genes <- as.integer(genes - 1)
  } else if (length(genes) == 1 && genes == "all_genes") {
    gene_names <- colnames(data)
  } else if (length(genes) == 1 && genes == "pca_only") {
    gene_names <- paste0("PC", 1:npca)
  } else {
    # character vector
    if (!all(genes %in% colnames(data))) {
      warning(paste0("Genes ", genes[!(genes %in% colnames(data))],
                     " not found.", collapse=", "))
    }
    genes <- which(colnames(data) %in% genes)
    gene_names <- colnames(data)[genes]
    genes <- as.integer(genes - 1)
  }

  # store parameters
  params <- list("data" = data, "k" = k, "alpha" = alpha, "t" = t,
                 "npca" = npca, "knn.dist.method" = knn.dist.method)
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
                          verbose = verbose)
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
                              verbose = verbose)
  }
  result <- operator$fit_transform(data,
                                   genes = genes,
                                   t_max = t.max)
  result <- as.data.frame(result)
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
#' data(magic_testdata)
#' data_magic <- magic(magic_testdata)
#' print(data_magic)
#' ## MAGIC with elements
#' ## $result : (500, 197)
#' ## $operator : Python MAGIC operator
#' ## $params : list with elements (data, k, alpha, t, npca, knn.dist.method)
#'
#' }
#' @rdname print
#' @method print magic
#' @export
print.magic <- function(x, ...) {
  result <- paste0("MAGIC with elements\n",
                   "  $result : (", nrow(x$result), ", ",
                   ncol(x$result), ")\n",
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
#' if (reticulate::py_module_available("magic")) {
#'
#' data(magic_testdata)
#' data_magic <- magic(magic_testdata)
#' summary(data_magic)
#' ## ZEB1
#' ## Min.   :0.01071
#' ## 1st Qu.:0.01119
#' ## Median :0.01130
#' ## Mean   :0.01129
#' ## 3rd Qu.:0.01140
#' ## Max.   :0.01201
#'
#' }
#' @rdname summary
#' @method summary magic
#' @export
summary.magic <- function(object, ...) {
  summary(object$result)
}

#' Convert a MAGIC object to a matrix
#'
#' Returns the smoothed data matrix
#'
#' @param x A fitted MAGIC object
#' @param ... Arguments for as.matrix()
#' @rdname as.matrix
#' @method as.matrix magic
#' @export
as.matrix.magic <- function(x, ...) {
  as.matrix(as.data.frame(x))
}
#' Convert a MAGIC object to a data.frame
#'
#' Returns the smoothed data matrix
#'
#' @param x A fitted MAGIC object
#' @param ... Arguments for as.data.frame()
#' @rdname as.data.frame
#' @method as.data.frame magic
#' @export
as.data.frame.magic <- function(x, ...) {
  x$result
}


#' Convert a MAGIC object to a data.frame for ggplot
#'
#' Passes the smoothed data matrix to ggplot
#' @importFrom ggplot2 ggplot
#' @param data A fitted MAGIC object
#' @param ... Arguments for ggplot()
#' @examples
#' if (reticulate::py_module_available("magic") && require(ggplot2)) {
#'
#' data(magic_testdata)
#' data_magic <- magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"))
#' ggplot(data_magic, aes(VIM, CDH1, colour=ZEB1)) +
#'   geom_point()
#'
#' }
#' @rdname ggplot
#' @method ggplot magic
#' @export
ggplot.magic <- function(data, ...) {
  ggplot2::ggplot(as.data.frame(data), ...)
}
