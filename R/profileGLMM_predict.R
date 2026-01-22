#' @title Prediction of cluster memberships and outcomes (internal function for predict method)
#' @description (This documentation is now for internal use only)
#' @param object An object of class \code{pglmm_fit} .
#' @param XFE A numeric matrix of fixed effects covariates for the prediction data.
#' @param XLat A numeric matrix of latent effect covariates.
#' @param UCont A numeric matrix or vector of continuous profile variables. Defaults to \code{NULL}.
#' @param UCat A numeric matrix or vector of categorical profile variables. Defaults to \code{NULL}.
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats dnorm
#' @importFrom Matrix KhatriRao t sparse.model.matrix
#' @noRd
#'
profileGLMM_predict = function(post_Obj, XFE, XLat, UCont, UCat){
  if (!inherits(post_Obj, "pglmm_fit")) {
    stop("object must be of class 'pglmm_fit'")
  }
  pred = {}
  n = dim(XFE)[1]
  pred$FE = as.matrix(XFE)%*%t(post_Obj$pop$betaFE['mean',])
  pred$Y = pred$FE

  if(!is.null(post_Obj$clust)){
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
    gamVec = c(post_Obj$clust$gamma[,sort(unique(pred$classPred))])
    pred$Int = t(KhatriRao(t(sparse.model.matrix(~ 0 + factor(pred$classPred))),t(XLat),make.dimnames = TRUE))%*%gamVec
    pred$Y = pred$Y + pred$Int
  }else{
    warning('No representative clustering provided')
    pred$classPred = NULL
    pred$Int = NULL
  }


  return(pred)

}




