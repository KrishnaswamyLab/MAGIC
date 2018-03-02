#' @title Do the magic!
#'
#' @param data the data in some format
#' @param t_diffusion diffusion time; suggested to start with 6 and increase it
#' up to 12. When no t_diffusion is given, a suggested optimal value of t is
#' calculated automatically
#' @param lib_size_norm Default: TRUE.
#' @param log_transform Default: FALSE.
#' @param pseudo-count A number indicating how much should be added to avoid
#'  log-transformation of zeros. Default: 0.1.
#' @param npca number of PCA components that should be used; default: 20.
#' @param k Number of nearest neighbors to use when running MAGIC. Default: 30.
#' @param ka kNN-autotune parameter for running MAGIC; default: 4
#' @param epsilon a value for the standard deviation of the kernel. Default: 1.
#' \code{epsilon = 0} is the uniform kernel.
#' @param rescale_percent To which percentile should the data be re-scaled.
#' Note: Do not set this higher than 0 if you also want to log-transform.
#' Default: 0.
#'
#' @export
#'
run_magic <- function(data, t_diffusion=0, lib_size_norm=TRUE,
                      log_transform=FALSE,
                      pseudo_count=0.1,
                      npca=100, k=12,
                      ka=4, epsilon=1, rescale_percent=0) {

  if (lib_size_norm){
    print('Library size normalization')
    libsize <- rowSums(data)
    data <- data / libsize * median(libsize)
  }

  if (log_transform){
    print('Log transforming')
    data <- log(data + pseudo_count)
  }

  N = nrow(data) # number of cells

  print('PCA')
  data_centr <- scale(data, scale=FALSE)
  svd <- rsvd::rsvd(t(data_centr), k=npca)
  data_pc <- data_centr %*% svd$u

  print('Computing distances')
  knndata <- FNN::get.knn(data_pc, k=k)
  idx <- knndata$nn.index
  dist <- knndata$nn.dist

  if (ka > 0) {
    print('Adapting sigma')
    dist <- dist / dist[,ka]
  }

  i <- rep((1:N), k)
  j <- c(idx)
  s <- c(dist)
  if (epsilon > 0) {
    W <- Matrix::sparseMatrix(i, j, x=s)
  } else {
    W <- Matrix::sparseMatrix(i, j, x=1) # unweighted kNN graph
  }

  print('Symmetrize distances')
  W <- as.matrix(W)
  W <- W + t(W)

  if (epsilon > 0){
    print('Computing kernel')
    W <- as(W, "sparseMatrix")
    Q <- Matrix::summary(W)
    i <- Q$i
    j <- Q$j
    s <- Q$x
    s <- s / epsilon^2
    s <- exp(-s);
    W <- Matrix::sparseMatrix(i, j, x=s)
  }

  print('Markov normalization')
  W <- as.matrix(W) # to dense matrix
  W <- W / rowSums(W) # Markov normalization

  print('Diffusing')
  if (t_diffusion == 0) {
    t_diffusion <- compute_optimal_t(data, W, t_max=12, n_genes=500)
  }
  W_t <- expm::"%^%"(W, t_diffusion)

  print('Imputing')
  data_imputed <- W_t %*% as.matrix(data)

  if (rescale_percent > 0) {
    print('Rescaling')

    M99 <- apply(data, 2, function(x) quantile(x, rescale_percent))
    M100 <- apply(data, 2, max)
    indices <- which(M99 == 0, arr.ind=TRUE)
    if (length(indices) > 0){
      M99[indices] <- M100[indices]
    }

    M99_new <- apply(data_imputed, 2, function(x) quantile(x, rescale_percent))
    M100_new <- apply(data_imputed, 2, max)
    indices <- which(M99_new == 0, arr.ind=TRUE)
    if (length(indices) > 0) {
      M99_new[indices] <- M100_new[indices]
    }

    max_ratio <- M99 / M99_new
    data_imputed <- data_imputed * matrix(rep(max_ratio, length(data[,1])),
                                          nrow=length(data[,1]), byrow=TRUE)
  }

  return(data.frame(data_imputed))

}

#' @description Compute optimal t automatically
#'
#' @param data input data
#' @param diff_op diffusion operator
#' @param t_max maximum number of t
#' @param n_genes number of genes
#' @param make_plots create a plot of R2 with respect to t
#' @return the optimal t and a vector of R2 values for all t
compute_optimal_t <- function(data, diff_op,
                             t_max=32, n_genes=ncol(data),
                             make_plots=TRUE) {
  idx_genes <- sample(1:ncol(data), n_genes)
  data_imputed <- data[,idx_genes]
  data_imputed <- data.matrix(data_imputed)

  if (min(data_imputed) < 0) {
    print('Data has negative values, shifting to positive')
    data_imputed = data_imputed - min(data_imputed)
  }

  r2_vec <- rep(NA, t_max)
  data_prev <- data_imputed / sum(data_imputed)
  print('Computing optimal t')
  for (i in 1:t_max) cat('.')
  cat('\r')
  for (i in 1:t_max) {
    data_imputed <- diff_op %*% data_imputed
    data_curr <- data_imputed / sum(data_imputed)
    cvec <- c(t(data_curr)) # Unroll the data frame to a vector.
    pvec <- c(t(data_prev))
    r2 <- r_square(cvec, pvec)[1]
    r2_vec[i] <- 1 - r2
    data_prev <- data_curr
    cat('*')
  }
  cat('\n')
  t_opt <- min(which(r2_vec < 0.05)) + 1
  print(paste('Optimal t =', t_opt))
  if (make_plots) {
    pdf("optimal_t.pdf", width=6, height=4)
    plot(r2_vec, type='l', xlab='t', ylab='R2')
    title('R2 vs. t')
    axis(side=1, at=1:length(r2_vec))
    abline(h=0.05, lty=2)
    points(t_opt, r2_vec[t_opt], col='red')
    dev.off()
  }
  return(t_opt)
}

#' @description Compute coefficient of determination of data fit model and RMSE
#'
#' @param y actual data
#' @param f model fit
#' @param c constant term in model
#' R-square may be a questionable measure of fit when no constant term is
#' included in the model.
#' [DEFAULT] TRUE : Use traditional R-square computation
#'          FALSE : Uses alternate R-square computation for model
#'                  without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
#' @return the coefficient of determination and the root mean squared error
r_square <- function(y, f, c=TRUE) {
  if (!all(dim(y) == dim(f))) {
    stop('Y and F must be the same size');
  }

  # Flatten array.
  y = c(y)
  f = c(f)

  # Remove NaN values.
  not_NaN <- !(is.nan(y) | is.nan(f))
  y <- y[not_NaN]
  f <- f[not_NaN]

  if (c) {
    r2 = max(0, 1 - sum((y - f)^2) / sum((y - mean(y))^2))
  } else {
    r2 = 1 - sum((y - f)^2) / sum((y)^2)
    if (r2 < 0) {
      warning('Consider adding a constant term to your model')
      r2 <- 0
    }
  }
  rmse = sqrt(mean((y - f)^2))
  return(c(r2, rmse))
}

