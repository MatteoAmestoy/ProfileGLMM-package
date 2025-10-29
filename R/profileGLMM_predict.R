#' Prediction of cluster memberships and post_Objcomes
#'
#' @param MCMC_Obj Profile GLMM  MCMC post_Objput of  profileGLMM_Gibbs function
#' @returns The Gibbs-sampled posterior. post_Objput of the cpp GSLoopCPP function.
#' @export
#'
#' @examples
profileGLMM_predict = function(post_Obj, XFE, XLat, UCont, UCat){
  pred= {}
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
    pred$Lat = t(KhatriRao(t(t(as(factor(pred$classPred),Class = "sparseMatrix"))),t(XLat),make.dimnames = TRUE))%*%gamVec
    pred$Y = pred$Y + pred$Lat
  }else{
    print('No representative clustering provided')
    pred$classPred = NULL
    pred$Lat = NULL
  }


  return(pred)

}




