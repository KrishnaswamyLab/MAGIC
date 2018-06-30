# Return TRUE if x and y are equal or both NA
null_equal <- function(x, y) {
  if (is.null(x) && is.null(y)) {
    return(TRUE)
  } else if (is.null(x) || is.null(y)) {
    return(FALSE)
  } else {
    return(x == y)
  }
}

load_pymagic <- function() {
  pymagic <<- reticulate::import("magic", delay_load = TRUE)
}

#' Install PHATE Python Package
#'
#' Install PHATE Python package into a virtualenv or conda env.
#'
#' On Linux and OS X the "virtualenv" method will be used by default
#' ("conda" will be used if virtualenv isn't available). On Windows,
#' the "conda" method is always used.
#' As of reticulate v1.7, this functionality is only available in the
#' development version of reticulate, which can be installed using
#' `devtools::install_github('rstudio/reticulate')`
#'
#' @param envname Name of environment to install packages into
#' @param method Installation method. By default, "auto" automatically finds
#' a method that will work in the local environment. Change the default to
#' force a specific installation method. Note that the "virtualenv" method
#' is not available on Windows.
#' @param conda Path to conda executable (or "auto" to find conda using the PATH
#'  and other conventional install locations).
#' @param pip Install from pip, if possible.
#' @param ... Additional arguments passed to conda_install() or
#' virtualenv_install().
#'
#' @export
install.magic <- function(envname = "r-reticulate", method = "auto",
                          conda = "auto", pip=TRUE, ...) {
  stop(paste0(
      "Cannot install MAGIC, please install from a console with ",
      "pip install --user git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python"
    ))
}

pymagic <- NULL

.onLoad <- function(libname, pkgname) {
  load_pymagic()
}
