
#' Initialise the priors of the model
#'
#' @param params A list with Fields: \itemize{
#' \item TBW}
#' @param nC An integer, the number of possible clusters
#'
#' @returns A list of the priors parameters for each variable
#'
#' @examples
prior_init = function(params, nC){
  ' initialise the priors of the model
-------Input:
      - params: list with Fields:
                    - TBW
      - priorInit (optionnal): list with fields
                    - TBW
-------Output:
  Note:
  '
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
