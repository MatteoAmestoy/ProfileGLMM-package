#' Initatilize the variables for the gibbs sampler.
#'
#' @param prior A list of prior configuration to draw initialization from.
#' Same structure as the output of prior_init.
#' @param params A list of the problem parameters.
#' Same structure as the output of process_Data_outcome.
#' @param nC number of possible clusters
#'
#' @returns A list containing the initalisation values of the Gibbs sampler
#' @export
#'
#' @examples
theta_init = function(prior,params,nC){

  theta = {}
  theta$sig2 = rgamma(1,1,1)#rgamma(1,prior$FE$a,prior$FE$b)

  theta$betaFE = rnorm(params$qFE,0,theta$sig2/sqrt(prior$FE$lambda))

  if(params$qRE>0){
    theta$SigRE = rinvwishart(prior$RE$eta,prior$RE$Phi) #
  }else{
    theta$SigRE = diag(1)
  }

  theta$SigLat = rinvwishart(prior$Lat$eta,prior$Lat$Phi) #
  theta$betaLat = rnorm(params$qLat,0,theta$sig2/sqrt(prior$FE$lambda))

  theta$gammaLat = matrix(0,nrow = params$qLat,ncol = nC)
  for (c in 1:nC){
    theta$gammaLat[,c] = rmvn(1,rep(0,params$qLat),theta$SigLat)}

  theta$SigmaClus = array(0, dim = c(params$qU,params$qU,nC))
  for (c in 1:nC){
    theta$SigmaClus[,,c] = rinvwishart(prior$assign$nu,prior$assign$Psi)}

  theta$muClus = matrix(0,nrow = params$qU,ncol = nC)
  for (c in 1:nC){
    theta$muClus[,c] = rmvn(1,prior$assign$mu,theta$SigmaClus[,,c]/prior$assign$lambda)*0}

  return(theta)
}
