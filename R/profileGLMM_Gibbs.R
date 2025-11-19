#' @title R Wrapper for Profile GLMM Gibbs Sampler (C++ backend)
#'
#' @description This is the main function for fitting the Profile Generalized Linear Mixed Model (Profile GLMM) using a blocked Gibbs sampling algorithm. It acts as an R wrapper, passing pre-processed data, initial values, and prior hyperparameters contained in the \code{model} object directly to the C++ implementation \code{GSLoopCPP}. The function simulates the posterior distribution of all model parameters, including fixed effects, random effects variance, profile cluster parameters, latent effects, and cluster assignments.
#'
#' @param model A list object containing all data, initial parameter values, model dimensions, prior hyperparameters, and model configuration (e.g., regression type). This object is typically the output of a data processing function like \code{process_Data_outcome}. Key components include:
#' \itemize{
#'   \item{\code{d}:}{ Data matrices (Y, XFE, XRE, XLat, UCont, UCat).}
#'   \item{\code{params}:}{ Model dimension parameters (e.g., nC, qRE, qUCont).}
#'   \item{\code{theta}:}{ Initial values for parameters ($\beta_{FE}, \sigma^2, \Sigma_{RE}, \text{cluster means}, \text{cluster covariance}, \text{cluster prob. vectors}, \Sigma_{Lat}, \gamma_{Lat}$).}
#'   \item{\code{prior}:}{ Hyperparameters for all prior distributions (e.g., $N, \text{Inverse-Wishart}, \text{Dirichlet}$).}
#'   \item{\code{regType}:}{ The type of regression being performed.}
#' }
#' @param nIt Integer, the total number of MCMC iterations *counting* the burn-in period. The sampler will run for \code{nIt - nBurnIn} iterations in total.
#' @param nBurnIn Integer, the number of initial MCMC iterations that are discarded (not saved) to allow the chain to converge.
#' @returns A list containing the saved Gibbs-sampled MCMC chains for all model parameters (e.g., \code{beta}, \code{Z}, \code{gamma}, \code{pvec}, \code{muClus}, \code{PhiClus}, etc.) and the variable names from the original data. This output is ready for post-processing with \code{profileGLMM_postProcess}.
#' @export
#'
#' @examples
#' # Assuming 'dataProfile' is the output of the data pre-processing step
#' MCMC_Obj = profileGLMM_Gibbs(model = dataProfile,
#' #   nIt = 5000,
#' #   nBurnIn = 1000
#' # )
profileGLMM_Gibbs = function(model,nIt,nBurnIn){


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
  return(gibbs_out)
}
