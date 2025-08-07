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

  gibbs_out = GSLoopCPP(nIt, nBurnIn,
                        proLMMObj$params$nC,
                        proLMMObj$d$Y,
                        proLMMObj$d$XFE,
                        proLMMObj$d$XRE,
                        proLMMObj$d$XLat,
                        proLMMObj$d$U,
                        proLMMObj$d$ZRE,
                        proLMMObj$theta$betaFE,
                        proLMMObj$theta$sig2,
                        proLMMObj$theta$SigRE,
                        proLMMObj$theta$muClus,
                        proLMMObj$theta$SigmaClus,
                        proLMMObj$theta$SigLat,
                        proLMMObj$theta$gammaLat,
                        proLMMObj$prior$FE$a,
                        proLMMObj$prior$FE$b,
                        proLMMObj$prior$FE$lambda,
                        proLMMObj$prior$RE$Phi,
                        proLMMObj$prior$RE$eta,
                        proLMMObj$prior$assign$lambda,
                        proLMMObj$prior$assign$mu,
                        proLMMObj$prior$assign$nu,
                        proLMMObj$prior$assign$Psi,
                        proLMMObj$prior$Lat$Phi,
                        proLMMObj$prior$Lat$eta,
                        proLMMObj$prior$DP$scale,
                        proLMMObj$prior$DP$shape)


  return(gibbs_out)
}
