#' Print method for pglmm_fit
#'
#' @param x An object of class \code{pglmm_fit}
#' @param ... Additional arguments
#'
#' @exportS3Method print pglmm_fit
print.pglmm_fit <- function(x, ...) {
  cat("--- pglmm_fit object containing the postprocessing of the MCMC samples ---\n")
  # Use a safe check (like is.null) to prevent errors during printing
  if (!is.null(x$clust$Kstar)) {
    cat("-- Representative clustering characteristics: \n")
    cat("Method used:", x$clust$mode, "\n")
    cat("Clusters found:", x$clust$Kstar, "\n")
  } else {
    cat("Representative clustering not computed.\n")
  }
  invisible(x)
}

#' @title Prediction of cluster memberships and outcomes
#' @description (This documentation is now for internal use only)
#' @param object An object of class \code{pglmm_fit} .
#' @param newData : A list with fields\itemize{
#' \item XFE A numeric matrix of fixed effects covariates for the prediction data.
#' \item XLat A numeric matrix of latent effect covariates.
#' \item UCont A numeric matrix or vector of continuous profile variables. Defaults to \code{NULL}.
#' \item UCat A numeric matrix or vector of categorical profile variables. Defaults to \code{NULL}.}
#' @param ... Additional arguments
#'
#' @exportS3Method predict pglmm_fit
#'
#' @examples
#' # Load MCMC_Obj, the result of profileGLMM_Gibbs()
#' data("examp")
#' MCMC_Obj = examp$MCMC_Obj
#'
#' # Post-process the results
#' post_Obj = profileGLMM_postProcess(MCMC_Obj, modeClus='LS')
#'
#'
predict.pglmm_fit <- function(x, newData, ...) {
  profileGLMM_predict(x, newData$XFE, newData$XLat, newData$UCont, newData$UCat)
}
