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

#' Print method for pglmm_fit
#'
#' @param x An object of class \code{pglmm_fit}
#' @param ... Additional arguments
#'
#' @exportS3Method print pglmm_fit
summary.pglmm_fit <- function(x, ...) {
  cat("--- pglmm_fit summary ---\n")

  cat("-- Fixed effect estimates: \n")
  print(x$pop$betaFE)
  cat(" \n")
  if (!is.null(x$clust$Kstar)) {
    cat("-- Representative clustering estimates -- \n")
    cat("Method used:", x$clust$mode,", clusters found:", x$clust$Kstar, "\n")
    cat("Cluster centroids: \n")
    print(x$clus$cen)
    cat("Cluster interaction parameters: \n")
    print(x$clus$gamma)
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
#' # Load post_Obj, the result of profileGLMM_postProcess()
#' data("examp")
#' post_Obj = examp$post_Obj
#'
#' # run prediction for training data
#' pred_Obj = predict(post_Obj,examp$dataProfile)
#'
#'
predict.pglmm_fit <- function(x, newData, ...) {
  return(profileGLMM_predict(x, newData$XFE, newData$XLat, newData$UCont, newData$UCat))
}
