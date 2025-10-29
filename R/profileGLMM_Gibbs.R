#' R wrapper to call the cpp GSLoopCPP function.
#'
#' @param model Profile LMM  preprocessed by the process_Data_outcome function
#' @param nIt Integer, total number of elements sampled by the gibbs sampler
#' @param nBurnIn Integer, number of elements of the chain that are not saved as part of the burn-in.
#'
#' @returns The Gibbs-sampled posterior. Output of the cpp GSLoopCPP function.
#' @export
#'
#' @examples
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
