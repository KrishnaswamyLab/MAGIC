#' @title Do the magic!
#'
#' @param data the data in some format
#' @param t_diffusion diffusion time; suggested to start with 6 and increase it
#' up to 12.
#' @param lib_size_norm Default: TRUE.
#' @param log_transform Default: FALSE.
#' @param pseudo-count A number indicating how much should be added to avoid
#'  log-transformation of zeros. Default: 0.1.
#' @param npca number of PCA components that should be used; default: 20.
#' @param k Number of nearest neighbors to use when running MAGIC. Default: 30.
#' @param ka kNN-autotune parameter for running MAGIC; default: 10
#' @param epsilon a value for the standard deviation of the kernel. Default: 1.
#' \code{epsilon = 0} is the uniform kernel.
#' @param rescale_percent To which percentile should the data be re-scaled.
#' Note: Do not set this higher than 0 if you also want to log-transform.
#' Default: 0.
#'
#' @export
#'
run_magic <- function(data, t_diffusion, lib_size_norm=TRUE,
                      log_transform=FALSE,
                      pseudo_count=0.1,
                      npca=20, k=30,
                      ka=10, epsilon=1, rescale_percent=0) {

  if (lib_size_norm){
    print('Library size normalization')
    libsize  = rowSums(data)
    data <- data / libsize * median(libsize)
  }

  if (log_transform){
    print('Log transforming')
    data <- log(data + pseudo_count)
  }

  N = nrow(data) # number of cells

  print('PCA')
  data_centr <- scale(data, scale = FALSE)
  svd <- rsvd(t(data_centr), k=npca)
  data_pc <- data_centr %*% svd$u

  print('Computing distances')
  knndata <- get.knn(data_pc, k=k)
  idx <- knndata$nn.index
  dist <- knndata$nn.dist

  if (ka>0)
    print('Adapting sigma')
    dist <- dist / dist[,ka]

  i <- rep((1:N), ncol(idx))
  j <- c(idx)
  s <- c(dist)
  if (epsilon > 0)
    W <- sparseMatrix(i, j, x = s)
  else
    W <- sparseMatrix(i, j, x = 1) # unweighted kNN graph

  print('Symmetrize distances')
  W <- as.matrix(W)
  W <- W + t(W);

  if (epsilon > 0){
    print('Computing kernel')
    Q <- summary(W)
    i <- Q$i
    j <- Q$j
    s <- Q$x
    s <- s / epsilon^2
    s <- exp(-s);
    W <- sparseMatrix(i, j, x = s)
  }

  print('Markov normalization')
  W <- W / rowSums(W) # Markov normalization

  W <- as.matrix(W) # to dense matrix

  print('Diffusing')
  W_t <- W^t

  print('Imputing')
  data_imputed <- W_t %*% as.matrix(data)

  if (rescale_percent > 0){
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
    data_imputed <- data_imputed * matrix(rep(max_ratio,length(data[,1])),nrow=length(data[,1]), byrow=TRUE)

  }
  return(data.frame(data_imputed))

}

#' Another function
#'
#' @description Not sure what it does.
#'
#' @param A some parameter
#' @param n some other parameter
#' @return not sure
"%^%"<-function(A,n){
  if(n==1) A else {B<-A; for(i in (2:n)){A<-A%*%B}}; A
}
