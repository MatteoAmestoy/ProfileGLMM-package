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
  # theta$betaLat = rnorm(params$qLat,0,theta$sig2/sqrt(prior$FE$lambda))

  theta$gammaLat = matrix(0,nrow = params$qLat,ncol = nC)
  for (c in 1:nC){
    theta$gammaLat[,c] = rmvn(1,rep(0,params$qLat),theta$SigLat)}


  theta$ClusCont = {}
  if(params$qUCont>0){
    theta$ClusCont$Sigma = array(0, dim = c(params$qUCont,params$qUCont,nC))
    for (c in 1:nC){
      theta$ClusCont$Sigma[,,c] = rinvwishart(prior$assign$Cont$nu,prior$assign$Cont$Psi)
    }
    theta$ClusCont$mu = matrix(0,nrow = params$qUCont,ncol = nC)
    for (c in 1:nC){
      theta$ClusCont$mu[,c] = rmvn(1,prior$assign$Cont$mu,theta$ClusCont$Sigma[,,c]/prior$assign$Cont$lambda)}
  } else{
    theta$ClusCont$Sigma = diag(1)
    theta$ClusCont$mu = matrix(0,nrow = 1,ncol = nC)
  }

  theta$ClusCat$pvecClus = matrix(0,nrow = length(params$catInd),ncol=nC)
  if(params$qUCat>0){
    tracker = 0
    for (cat in 1:params$qUCat){
      n_ = sum(params$catInd==(cat-1))
      theta$ClusCat$pvecClus[(1:n_)+tracker,]= rdirichlet(nC, rep(prior$assign$Cat$alpha[cat],n_))
      tracker = tracker + n_
    }

  }

  return(theta)
}
