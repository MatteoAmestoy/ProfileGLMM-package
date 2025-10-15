#' R wrapper to call the cpp GSLoopCPP function.
#'
#' @param proLMMObj Profile LMM  preprocessed by the process_Data_outcome function
#' @param nIt Integer, total number of elements sampled by the gibbs sampler
#' @param nBurnIn Integer, number of elements of the chain that are not saved as part of the burn-in.
#'
#' @returns The Gibbs-sampled posterior. Output of the cpp GSLoopCPP function.
#' @export
#'
#' @examples
profileGLMM_Gibbs = function(proLMMObj,nIt,nBurnIn){


  if(proLMMObj$regType==1){proLMMObj$theta$sig2=1}
  if (is.null(proLMMObj$d$XRE)){proLMMObj$d$XRE = matrix() }

  gibbs_out = GSLoopCPP(nIt, nBurnIn,
                        proLMMObj$params$nC,
                        proLMMObj$params$qRE,
                        proLMMObj$params$qUCont, #-------------
                        proLMMObj$d$Y,
                        proLMMObj$d$XFE,
                        proLMMObj$d$XRE,
                        proLMMObj$d$XLat,
                        proLMMObj$d$UCont,
                        proLMMObj$d$UCat,
                        proLMMObj$params$catInd,
                        proLMMObj$d$ZRE,
                        proLMMObj$theta$betaFE,
                        proLMMObj$theta$sig2,
                        proLMMObj$theta$SigRE,
                        proLMMObj$theta$ClusCont$mu,
                        proLMMObj$theta$ClusCont$Sigma,
                        proLMMObj$theta$ClusCat$pvecClus,#mat
                        proLMMObj$theta$SigLat,
                        proLMMObj$theta$gammaLat,
                        proLMMObj$prior$FE$a,
                        proLMMObj$prior$FE$b,
                        proLMMObj$prior$FE$lambda,
                        proLMMObj$prior$RE$Phi,
                        proLMMObj$prior$RE$eta,
                        proLMMObj$prior$assign$Cont$lambda,
                        proLMMObj$prior$assign$Cont$mu,
                        proLMMObj$prior$assign$Cont$nu,
                        proLMMObj$prior$assign$Cont$Psi,
                        proLMMObj$prior$assign$Cat$alpha,#vec
                        proLMMObj$prior$Lat$Phi,
                        proLMMObj$prior$Lat$eta,
                        proLMMObj$prior$DP$scale,
                        proLMMObj$prior$DP$shape,
                        proLMMObj$regType)

  return(gibbs_out)
}
