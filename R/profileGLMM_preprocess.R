
#' Preprocess the data from a list describing the profile LMM model
#'
#' @param covList A list with fields:\itemize{
#' \item  FE fixed effect covariates names/index in dataframe
#' \item  RE random effect covariates names/index in dataframe
#' \item  Lat latent effect covariates names/index in dataframe
#' \item  Assign assignement variable categorical not supported yet
#' \item  REunit statistical unit of the RE colomn name/index
#' \item  Y outcome (Continuous)}
#' @param dataframe A dataframe containing outcome anf covariates
#' @param nC int: maximal number of cluster for the DP truncation
#' @param intercept (optionnal): A list with fields\itemize{
#' \item  RE bool indicating if FE have an intercept
#' \item  FE bool indicating if RE have an intercept
#' \item  Lat bool indicating if Latent have an intercept}
#'
#' @returns A list with\itemize{
#' \item  d dictionary with  [XFE,XRE,XLat,U,ZRE] design matrices
#' \item  [[params]] list of the parameters of the data\itemize{
#'                  \item  n int nb of obs
#'                  \item  qFE lint, number of covariates of FE
#'                  \item  nRE int, number of stat units of RE
#'                  \item  qRE int, number of covariates of RE}
#' \item prior a list with all the specification of the default prior used
#' \item theta a list with a default set of parameters to start the chain, drawn from the prior}
#' @export
#'
#' @examples
profileGLMM_preprocess <- function(covList, dataframe, nC, intercept = list(FE=T,RE=T,Lat =T)) {
  ' Preprocess the data from a LMER formula
-------Input:
      - dataframe: dataframe on which to apply the formula
      - covList: list with Fields:
                    - FE fixed effect covariates names/index in dataframe
                    - RE random effect covariates names/index in dataframe
                    - Lat latent effect covariates names/index in dataframe
                    - Assign assignement variable categorical not supported yet
                    - REunit statistical unit of the RE colomn name/index
                    - Y outcome (Continuous)
      - intercept (optionnal): list with fields
                    - RE bool indicating if FE have an intercept
                    - FE bool indicating if RE have an intercept
                    - Lat bool indicating if Latent have an intercept
-------Output:
      - d dictionary with [XFE,XRE,XLat,U,ZRE] design matrices
      - [[params]] list of the parameters of the data
            - n int nb of obs
            - qFE lint, number of covariates of FE
            - nRE int, number of stat units of RE
            - qRE int, number of covariates of RE
  Note:
  '

  d = {}
  d$names = {}
  d$XFE = encodeCat(dataframe[,covList$FE, drop = FALSE])
  d$names$FE = colnames(d$XFE)
  if (intercept$FE){d$XFE = cbind(1,d$XFE)
  d$names$FE= c('Intercept',d$names$FE) }
  d$XFE = as.matrix(d$XFE)

  d$XRE = encodeCat(dataframe[,covList$RE, drop = FALSE])
  d$names$RE = colnames(d$XRE)
  if (intercept$RE){d$XRE = cbind(1,d$XRE)
  d$names$RE= c('Intercept',d$names$RE) }
  d$XRE = as.matrix(d$XRE)

  d$XLat = encodeCat(dataframe[,covList$Lat, drop = FALSE])
  d$names$Lat = colnames(d$XLat)
  if (intercept$Lat){d$XLat = cbind(1,d$XLat)
  d$names$Lat= c('Intercept',d$names$Lat) }
  d$XLat = as.matrix(d$XLat)

  d$U = dataframe[,covList$Assign, drop = FALSE]
  d$U = as.matrix(d$U)
  d$names$U = colnames(d$U)

  d$Y = drop(dataframe[,covList$Y])

  d$ZRE = as.numeric(factor(dataframe[,covList$REunit]))-1

  params = {}
  params$n = dim(d$XFE)[1]
  params$nRE = max(d$ZRE)+1
  params$qFE = dim(d$XFE)[2]
  params$qRE = dim(d$XRE)[2]
  params$qLat = dim(d$XLat)[2]
  params$qU = dim(d$U)[2]
  params$nC = nC

  prior = prior_init(params, nC)

  theta = theta_init(prior,params,nC)

  return(list(d = d,
              params = params,
              prior = prior,
              theta = theta))

}
