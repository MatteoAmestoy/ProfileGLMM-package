#' @title R Wrapper for Profile GLMM Gibbs Sampler (C++ backend)
#'
#' @description This is the main function for fitting the Profile Generalized Linear Mixed Model
#' using a blocked Gibbs sampling algorithm. It acts as an R wrapper, passing an object of class
#' \code{pglmm_data} directly to the RCPP implementation \code{GSLoopCPP}. The function simulates the
#' posterior distribution of all model parameters, including fixed effects, random effects variance,
#' profile cluster parameters, latent effects, and cluster assignments.
#'
#' @param model An object of class \code{glmm_data} (the output of \code{profileGLMM_preprocess}).
#' This contains the design matrices, initial values, dimensions, and prior hyperparameters.
#' @param nIt Integer, the total number of MCMC iterations counting the burn-in period.
#' The sampler will return \code{nIt - nBurnIn} iterations in total.
#' @param nBurnIn Integer, the number of initial MCMC iterations that are discarded (not saved)
#' to allow the chain to converge.
#'
#' @return An object of class \code{pglmm_mcmc}. This is a list containing the saved Gibbs-sampled
#' MCMC chains for all model parameters (e.g., \code{beta}, \code{Z}, \code{gamma}, \code{pvec}, \code{muClus}, \code{PhiClus}, etc.)
#' and the variable names from the original data.
#' This output is intended for post-processing with \code{profileGLMM_postProcess}.
#'
#' @export
#'
#' @examples
#' # Load examp, which contains a pre-processed pglmm_data object
#' data("examp")
#' dataProfile = examp$dataProfile
#'
#' # Run the Gibbs Sampler
#' MCMC_Obj = profileGLMM_Gibbs(
#'   model = dataProfile,
#'   nIt = 100,
#'   nBurnIn = 10
#' )
profileGLMM_Gibbs = function(model,nIt,nBurnIn){

  if (!inherits(model, "pglmm_data")) stop("Input must be a pglmm_data object.")

  if(model$regType==1){model$theta$sig2=1}
  if (is.null(model$d$XRE)){model$d$XRE = matrix() }
  if (model$params$qUCat == 0){model$d$UCat = matrix() }

  gibbs_out = GSLoopCPP(nIt+1, nBurnIn,
                        model$params$nC,
                        model$params$qRE,
                        model$params$qUCont, #-------------
                        model$d$Y,
                        model$d$XFE,
                        model$d$XRE,
                        model$d$XLat,
                        model$d$UCont,
                        model$d$UCat,
                        model$params$catInd,
                        model$d$ZRE,
                        model$theta$betaFE,
                        model$theta$sig2,
                        model$theta$SigRE,
                        model$theta$ClusCont$mu,
                        model$theta$ClusCont$Sigma,
                        model$theta$ClusCat$pvecClus,#mat
                        model$theta$SigLat,
                        model$theta$gammaLat,
                        model$prior$FE$a,
                        model$prior$FE$b,
                        model$prior$FE$lambda,
                        model$prior$RE$Phi,
                        model$prior$RE$eta,
                        model$prior$assign$Cont$lambda,
                        model$prior$assign$Cont$mu,
                        model$prior$assign$Cont$nu,
                        model$prior$assign$Cont$Psi,
                        model$prior$assign$Cat$alpha,#vec
                        model$prior$Lat$Phi,
                        model$prior$Lat$eta,
                        model$prior$DP$scale,
                        model$prior$DP$shape,
                        model$regType)
  gibbs_out$names = model$d$names
  class(gibbs_out) <- "pglmm_mcmc"
  return(gibbs_out)
}
