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

load_pymagic <- function(delay_load = FALSE) {
    if (is.null(pymagic)) {
    result <- try(pymagic <<- reticulate::import("magic", delay_load = delay_load))
  } else {
    result <- try(reticulate::import("magic", delay_load = delay_load))
  }
  if (methods::is(result, "try-error") &&
      (length(grep("ModuleNotFoundError: No module named 'magic'", result)) > 0 ||
        length(grep("ImportError: No module named magic", result)) > 0)) {
      if (utils::menu(c("Yes", "No"), title="Install MAGIC Python package with reticulate?") == 1) {
        install.magic()
      }
  } else if (length(grep("r\\-reticulate", reticulate::py_config()$python)) > 0) {
    message("Consider removing the 'r-reticulate' environment by running:")
    if (grep("virtualenvs", reticulate::py_config()$python)) {
      message("reticulate::virtualenv_remove('r-reticulate')")
    } else {
      message("reticulate::conda_remove('r-reticulate')")
    }
  }
}

#' Install MAGIC Python Package
#'
#' Install MAGIC Python package into a virtualenv or conda env.
#'
#' On Linux and OS X the "virtualenv" method will be used by default
#' ("conda" will be used if virtualenv isn't available). On Windows,
#' the "conda" method is always used.
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
  message("Attempting to install MAGIC python package with reticulate")
  tryCatch({
    reticulate::py_install("magic-impute",
      envname = envname, method = method,
      conda = conda, pip=pip, ...
    )
    message("Install complete. Please restart R and try again.")
  },
  error = function(e) {
    stop(paste0(
      "Cannot locate MAGIC Python package, please install through pip ",
      "(e.g. pip install magic-impute) and then restart R."
    ))
  }
  )
}

pymagic <- NULL

.onLoad <- function(libname, pkgname) {
  py_config <- reticulate::py_discover_config(required_module = "magic")
  load_pymagic(delay_load = TRUE)
}
