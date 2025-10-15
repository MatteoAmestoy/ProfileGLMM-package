
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
  prior$assign$Cont = {}
  prior$assign$Cont$lambda = 1
  prior$assign$Cont$mu = rep(0,max(params$qUCont,1))
  prior$assign$Cont$nu = max(params$qUCont,1)+1
  prior$assign$Cont$Psi = (max(params$qUCont,1)+prior$assign$Cont$nu+1)*diag(max(params$qUCont,1))
  prior$assign$Cat = {}
  prior$assign$Cat$alpha = rep(1,min(params$qUCat,1))


  prior$Lat = {}
  prior$Lat$eta = params$qLat+4
  prior$Lat$Phi = diag(params$qLat)*(prior$Lat$eta)/2


  prior$DP = {}
  prior$DP$scale = sqrt(nC)
  prior$DP$shape = sqrt(nC)


  return(prior)
}
