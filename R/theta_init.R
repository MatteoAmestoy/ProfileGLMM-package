#' @title Initialize the variables for the Gibbs sampler chain
#'
#' @description This function generates initial values (\code{theta}) for all parameters in the Profile GLMM Gibbs sampler by drawing from the specified prior distributions. These initial values are crucial for starting the MCMC chain in \code{profileGLMM_Gibbs}. The initialization includes parameters for fixed effects, random effects variance, latent effects, and the profile cluster parameters (centroids, covariances, and categorical probability vectors).
#'
#' @param prior A list containing the prior configuration to draw initialization from. This list should match the structure produced by the \code{prior_init} function, including hyperparameters for FE, RE, Latent, and cluster assignment priors.
#' @param params A list containing the problem's dimensional parameters and indices (e.g., number of observations, number of covariates). This list should match the structure of the output from \code{process_Data_outcome}.
#'
#' @returns A list (\code{theta}) containing the sampled initialization values for the Gibbs sampler. Key elements include:
#' \describe{
#'   \item{\code{sig2}:}{ Initial residual variance.}
#'   \item{\code{betaFE}:}{ Initial fixed effects coefficients.}
#'   \item{\code{SigRE}:}{ Initial random effects covariance matrix.}
#'   \item{\code{SigLat}:}{ Initial latent effects covariance matrix.}
#'   \item{\code{gammaLat}:}{ Initial latent effects coefficients, organized by cluster.}
#'   \item{\code{ClusCont}:}{ List containing initial continuous cluster parameters (\code{mu} and \code{Sigma}).}
#'   \item{\code{ClusCat}:}{ List containing initial categorical cluster parameters (\code{pvecClus}).}
#' }
#' @export
#' @importFrom LaplacesDemon rinvwishart rmvn
#' @importFrom stats rgamma rnorm
#' @importFrom MCMCpack rdirichlet
#'
#' @examples
#' \dontrun{
#' # Assuming prior_config and problem_params are available
#' # initial_theta <- theta_init(
#' #   prior = prior_config,
#' #   params = problem_params
#' # )
#' }
theta_init = function(prior,params){
  nC = params$nC
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
