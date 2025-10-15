
#' Preprocess the data from a list describing the profile LMM model
#' @param regType A string, current possibilities: linear or probit
#' @param covList A list with fields:\itemize{
#' \item  FE fixed effect covariates names/index in dataframe
#' \item  RE random effect covariates names/index in dataframe
#' \item  Lat latent effect covariates names/index in dataframe
#' \item  Assign assignement variables list with fields:\itemize{
#'        \item  Cont Continuous variables names/index in dataframe
#'        \item  Cat Categorical variables names/index in dataframe
#' }
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
#' \item  d dictionary with  [XFE,XRE,XLat,UCont,UCat,ZRE] design matrices
#' \item  [[params]] list of the parameters of the data\itemize{
#'                  \item  n int nb of obs
#'                  \item  qFE lint, number of covariates of FE
#'                  \item  nRE int, number of stat units of RE
#'                  \item  qRE int, number of covariates of RE
#'                  \item  qLat int, number of covariates interacting with the
#'                                   latent clusters
#'                  \item  qUCont int, number of continuous clustering covariates
#'                  \item  qUCat int, number of categorical clustering covariates
#'                  \item  nC int, maximal number of clusters}
#' \item prior a list with all the specification of the default prior used
#' \item theta a list with a default set of parameters to start the chain, drawn from the prior
#' \item regType an int. Currently 0 for linear, 1 for probit}
#' @export
#'
#' @examples
profileGLMM_preprocess <- function(regtype, covList, dataframe, nC, intercept = list(FE=T,RE=T,Lat =T)) {

  d = {}
  d$names = {}
  params = {}

  if (regtype == 'linear'){
    rT = 0
    d$Y = drop(dataframe[,covList$Y])
  }else if(regtype == 'probit'){
    rT = 1
    d$Y = factor(dataframe[,covList$Y])
    if( length(levels(d$Y))==2){
      print(paste0('Reference category (Y = 1) for Y = ',levels(d$Y)[2]))
      d$Y = (d$Y==levels(d$Y)[2])}else{
        stop(paste0(covList$Y,' has ',length(levels(d$Y)),' levels, while 2 expected.'))
      }
  }else{stop('Regression type not supported please choose out of linear or probit.')}

  n = length(d$Y)

  if (length(covList$FE)!=0){
    d$XFE = encodeCat(dataframe[,covList$FE, drop = FALSE])
    d$names$FE = colnames(d$XFE)
    if (intercept$FE){
      d$XFE = cbind(1,d$XFE)
      d$names$FE = c('Intercept',d$names$FE)}
    d$XFE = as.matrix(d$XFE)}else if(intercept$FE){
      d$names$FE = c('Intercept')
      d$XFE = as.matrix(rep(1,n))
    }else{
      stop('Error: no fixed effects provided')
    }


  if (length(covList$RE)!=0){
    d$XRE = encodeCat(dataframe[,covList$RE, drop = FALSE])
    d$names$RE = colnames(d$XRE)
    if (intercept$RE){d$XRE = cbind(1,d$XRE)
    d$names$RE= c('Intercept',d$names$RE) }
    d$XRE = as.matrix(d$XRE)}else if(intercept$RE){
      d$names$RE = c('Intercept')
      d$XRE = as.matrix(rep(1,n))
    }else{
      print('Warning: no random effects provided')
      d$XRE = NULL
      d$names$RE = NULL
    }

  if (length(covList$Lat)!=0){
    d$XLat = encodeCat(dataframe[,covList$Lat, drop = FALSE])
    d$names$Lat = colnames(d$XLat)
    if (intercept$Lat){d$XLat = cbind(1,d$XLat)
    d$names$Lat= c('Intercept',d$names$Lat) }
    d$XLat = as.matrix(d$XLat)}else if(intercept$Lat){
      d$names$Lat = c('Intercept')
      d$XLat = as.matrix(rep(1,n))
    }else{
      stop('ERROR: no cluster effect provided')
    }
  U_Empty = TRUE
  d$UCat = data.frame()
  if(length(covList$Assign$Cat)!=0){
    d$names$UCat = covList$Assign$Cat
    U_Empty = FALSE
    params$catInd = c()
    tracker = 1
    for(cat in covList$Assign$Cat){
      tmp = encodeCat(as,factor(dataframe[,cat, drop = FALSE]))
      d$UCat[,0:(dim(tmp)[2])+tracker] = cbind(1-rowSums(tmp,dims = 2),tmp)
      params$catInd = c(params$catInd,rep(cat,dim(d$UCat$cat)[2]))
      tracker = tracker + dim(d$UCat$cat)[2]
    }
  } else{
    params$catInd = c(-1)
  }
  if(length(covList$Assign$Cont)!=0){
    d$UCont = dataframe[,covList$Assign$Cont, drop = FALSE]
    d$names$UCont = colnames(d$UCont)
    d$UCont = as.matrix(d$UCont)
    U_Empty = FALSE
  }
  if(U_Empty){
    stop('ERROR: no clustering variables provided')}
  d$UCat = as.matrix(d$UCat)

  if (length(covList$REunit)!=0){
    d$ZRE = as.numeric(factor(dataframe[,covList$REunit]))-1
  }else {
    d$ZRE = as.matrix(rep(-1,n))
  }

  params$n = dim(d$XFE)[1]
  params$nRE = max(d$ZRE)+1
  params$qFE = dim(d$XFE)[2]
  if(!is.null(dim(d$XRE)[2])){
    params$qRE = dim(d$XRE)[2]} else{
      params$qRE = 0
    }
  if(!is.null(dim(d$UCont)[2])){
    params$qUCont = dim(d$UCont)[2]} else{
      params$qUCont = 0
    }
  if(dim(d$UCat)[2]!=0){
    params$qUCat = length(unique(params$catInd))} else{
      params$qUCat = 0
    }
  params$qLat = dim(d$XLat)[2]
  params$nC = nC

  prior = prior_init(params, nC)

  theta = theta_init(prior,params,nC)

  return(list(d = d,
              params = params,
              prior = prior,
              theta = theta,
              regType = rT))

}
