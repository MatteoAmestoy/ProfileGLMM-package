#' @title Initialize the prior hyperparameters for the Profile GLMM
#'
#' @description This function establishes the prior distributions for all parameters
#' in the Profile GLMM. It sets up vague, non-informative priors (often using small
#' precision/large variance or conjugate forms like Wishart/Dirichlet) for the fixed effects (\eqn{beta_{FE}}),
#' residual variance (\eqn{\sigma^2}), random effects covariance (\eqn{\Sigma_{RE}}), latent effects covariance (\eqn{\Sigma_{Lat}}),
#' cluster parameters (means and covariances), and the Dirichlet Process parameters (\eqn{\alpha}).
#'
#' @param params A list containing dimensional parameters of the model (often the output of \code{process_Data_outcome}). Important fields used for prior setup include:
#' \describe{
#'   \item{\code{qFE}:}{ Number of fixed effects coefficients.}
#'   \item{\code{qRE}:}{ Dimension of the random effects vector.}
#'   \item{\code{qLat}:}{ Dimension of the latent effects vector.}
#'   \item{\code{qUCont}:}{ Number of continuous profile variables.}
#'   \item{\code{qUCat}:}{ Number of categorical profile variables.}
#' }
#'
#' @returns A list (\code{prior}) containing the hyperparameter values structured by the parameter block they govern:
#' \describe{
#'   \item{\code{FE}:}{ Priors for fixed effects and residual variance (e.g., \code{lambda}, \code{a}, \code{b} for conjugate Normal-Gamma).}
#'   \item{\code{RE}:}{ Inverse-Wishart priors for random effects covariance (\eqn{\Sigma_{RE}}) (e.g., \code{Phi}, \code{eta}).}
#'   \item{\code{assign}:}{ Priors for the cluster assignment parameters, nested under \code{Cont} (Normal-Inverse-Wishart for continuous) and \code{Cat} (Dirichlet for categorical).}
#'   \item{\code{Lat}:}{ Inverse-Wishart prior for the latent effects covariance (\eqn{\Sigma_{Lat}}) (e.g., \code{Phi}, \code{eta}).}
#'   \item{\code{DP}:}{ Parameters for the Dirichlet Process prior (e.g., \code{scale}, \code{shape}).}
#' }
#' @export
#'
#' @examples
#' # Load dataProfile, the result of profileGLMM_preProcess()
#' data("examp")
#' dataProfile = examp$dataProfile
#' prior_config <- prior_init(dataProfile$params)

prior_init = function(params){
  nC = params$nC
  prior ={}

  prior$FE = {}
  prior$FE$lambda = 10**(-6)
  prior$FE$a = 10**(-6)
  prior$FE$b = 10**(-6)

  if(params$qRE>0){
    prior$RE = {}
    prior$RE$Phi = diag(params$qRE)
    prior$RE$eta = params$qRE}else{
      prior$RE = {}
      prior$RE$Phi = diag(1)
      prior$RE$eta = -1
    }

  prior$assign = {}
  prior$assign$Cont = {}
  prior$assign$Cont$lambda = 1
  prior$assign$Cont$mu = rep(0,max(params$qUCont,1))
  prior$assign$Cont$nu = max(params$qUCont,1)+1
  prior$assign$Cont$Psi = (max(params$qUCont,1)+prior$assign$Cont$nu+1)*diag(max(params$qUCont,1))
  prior$assign$Cat = {}
  prior$assign$Cat$alpha = rep(1,min(params$qUCat,1))


  prior$Lat = {}
  prior$Lat$eta = params$qLat
  prior$Lat$Phi = diag(params$qLat)


  prior$DP = {}
  prior$DP$scale = sqrt(nC)
  prior$DP$shape = sqrt(nC)


  return(prior)
}
