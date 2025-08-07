
#' Initialise the priors of the model
#'
#' @param params A list with Fields: \itemize{
#' \item FE fixed effect covariates names/index in dataframe
#' \item RE random effect covariates names/index in dataframe
#' \item REunit statistical unit of the RE colomn name/index
#' \item Lat latent effect covariates names/index in dataframe
#' \item Assign assignement variable categorical not supported yet}
#' @param nC An integer, the number of possible clusters
#'
#' @returns A list of the priors parameters for each variable
#'
#' @examples
prior_init = function(params, nC){
  ' initialise the priors of the model
-------Input:
      - params: list with Fields:
                    - FE fixed effect covariates names/index in dataframe
                    - RE random effect covariates names/index in dataframe
                    - REunit statistical unit of the RE colomn name/index
                    - Lat latent effect covariates names/index in dataframe
                    - Assign assignement variable categorical not supported yet
      - priorInit (optionnal): list with fields
                    - RE bool indicating if FE have an intercept
                    - FE bool indicating if RE have an intercept
                    - Lat bool indicating if Latent have an intercept
-------Output:
  Note:
  '
  prior ={}

  prior$FE = {}
  prior$FE$lambda = 10**(-6)
  prior$FE$a = 10**(-6)
  prior$FE$b = 10**(-6)

  prior$RE = {}
  prior$RE$Phi = diag(params$qRE)
  prior$RE$eta = params$qRE

  prior$assign = {}
  prior$assign$lambda = 1
  prior$assign$mu = rep(0,params$qU)
  prior$assign$nu = params$qU+1
  prior$assign$Psi = (params$qU+prior$assign$nu+1)*diag(params$qU)

  prior$Lat = {}
  prior$Lat$eta = params$qLat+4
  prior$Lat$Phi = diag(params$qLat)*(prior$Lat$eta)/2


  prior$DP = {}
  prior$DP$scale = sqrt(nC)
  prior$DP$shape = sqrt(nC)


  return(prior)
}
