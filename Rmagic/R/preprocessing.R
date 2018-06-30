#' Performs L1 normalization on input data such that the sum of expression
#' values for each cell sums to 1, then returns normalized matrix to the metric
#' space using median UMI count per cell effectively scaling all cells as if
#' they were sampled evenly.

#' @param data matrix (n_samples, n_dimensions)
#' 2 dimensional input data array with n cells and p dimensions
#' @param verbose boolean, default=FALSE. If true, print verbose output

#' @return data_norm matrix (n_samples, n_dimensions)
#' 2 dimensional array with normalized gene expression values
#' @import Matrix
#'
#' @export
library.size.normalize <- function(data, verbose=FALSE) {
  if (verbose) {
    message(paste0(
      "Normalizing library sizes for ",
      nrow(data), " cells"
    ))
  }
  library_size <- Matrix::rowSums(data)
  median_transcript_count <- stats::median(library_size)
  data_norm <- median_transcript_count * data / library_size
  data_norm
}
