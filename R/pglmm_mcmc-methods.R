#' Print method for pglmm_mcmc
#'
#' @param x An object of class \code{pglmm_mcmc}
#' @param ... Additional arguments
#'
#' @exportS3Method print pglmm_mcmc
print.pglmm_mcmc <- function(x, ...) {
  cat("--- pglmm_mcmc object containing", length(x$sig2),  " MCMC samples")
  invisible(x)
}
