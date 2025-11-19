#' @title Prediction of cluster memberships and outcomes
#'
#' @description This function uses the results of the post-processed Profile GLMM MCMC chain to predict cluster memberships and outcomes for new or existing data. It first calculates the fixed effect (FE) contribution and then, if a representative clustering is available in \code{post_Obj}, computes the predicted cluster membership and the corresponding latent effect (Lat) contribution to the outcome.
#'
#' @param post_Obj The post-processed output from the \code{profileGLMM_postProcess} function. Must contain \code{pop} for population constant parameters and optionally \code{clust} for cluster-specific parameters.
#' @param XFE A numeric matrix of fixed effects covariates for the prediction data.
#' @param XLat A numeric matrix of latent effect covariates. This matrix is used for the interaction term with the predicted cluster membership.
#' @param UCont A numeric matrix or vector of continuous profile variables (used for predicting cluster membership). Set to \code{NULL} if no continuous variables were used in the model.
#' @param UCat A numeric matrix or vector of categorical profile variables (used for predicting cluster membership). Set to \code{NULL} if no categorical variables were used in the model.
#' @returns A list with the following elements:
#' \itemize{
#'   \item{\code{FE}:}{ A numeric vector of the predicted fixed effects contribution to the outcome.}
#'   \item{\code{Y}:}{ A numeric vector of the total predicted outcome ($\text{FE} + \text{Lat}$).}
#'   \item{\code{classPred}:}{ A factor vector of the predicted cluster membership for each observation. \code{NULL} if no representative clustering was provided in \code{post_Obj}.}
#'   \item{\code{Int}:}{ A numeric vector of the predicted latent effect contribution to the outcome. \code{NULL} if no representative clustering was provided.}
#' }
#' @export
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats dnorm
#' @importFrom Matrix KhatriRao t
#'
#' @examples
#' \dontrun{
#' # Assuming post_Obj is the result of profileGLMM_postProcess()
#' pred_Obj = profileGLMM_predict(post_Obj,
#'                                dataProfile$d$XFE,
#'                                dataProfile$d$XLat,
#'                                dataProfile$d$UCont,
#'                                dataProfile$d$UCat)
#' }
profileGLMM_predict = function(post_Obj, XFE, XLat, UCont, UCat){
  pred = {}
  n = dim(XFE)[1]
  pred$FE = as.matrix(XFE)%*%t(post_Obj$pop$betaFE['mean',])
  pred$Y = pred$FE

  if(!is.null(post_Obj$clust)){
    gamVec = c(post_Obj$clust$gamma)
    matClassPred = matrix(1,nrow = n,ncol = post_Obj$clust$Kstar)
    for( c in 1:post_Obj$clust$Kstar){
      if(!is.null(UCont)){
        if (dim(UCont)[2]>1){
          matClassPred[,c] = matClassPred[,c] * dmvnorm(UCont,
                                                        post_Obj$clust$cen[,c],
                                                        drop(post_Obj$clust$coVar[,,c]),
                                                        log = F)
        }else{
          matClassPred[,c] = matClassPred[,c] * dnorm(UCont,
                                                      post_Obj$clust$cen[,c],
                                                      drop(post_Obj$clust$coVar[,,c]),
                                                      log = F)}

      }
      if(!is.null(UCat)){
        matClassPred[,c] = matClassPred[,c] * post_Obj$clust$pvec[c]}
    }

    pred$classPred = as.factor(apply(matClassPred,1,which.max))
    pred$Int = t(KhatriRao(as(factor(pred$classPred),Class = "sparseMatrix"),t(XLat),make.dimnames = TRUE))%*%gamVec
    pred$Y = pred$Y + pred$Int
  }else{
    print('No representative clustering provided')
    pred$classPred = NULL
    pred$Int = NULL
  }


  return(pred)

}




