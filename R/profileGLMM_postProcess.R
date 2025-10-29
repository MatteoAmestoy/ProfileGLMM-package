#' Post processing of the MCMC chain.
#'
#' @param MCMC_Obj Profile GLMM  MCMC output of  profileGLMM_Gibbs function
#' @returns The Gibbs-sampled posterior. Output of the cpp GSLoopCPP function.
#' @export
#'
#' @examples
profileGLMM_postProcess = function(MCMC_Obj, modeClus='NG', comp_cooc = T, alpha = 0.05){

  nSim = dim(MCMC_Obj$Z)[2]
  nUCat = dim(MCMC_Obj$pvec)[1]
  nUCont = dim(MCMC_Obj$muClus)[1]
  qLat = dim(MCMC_Obj$gamma)[1]


  # post processing of the population parameters
  pop = {}


  betaFE = apply(MCMC_Obj$beta, 1, mean)
  betaFECI = apply(MCMC_Obj$beta, 1, quantile,c(alpha/2,1-alpha/2))
  pop$betaFE = as.data.frame(rbind(betaFE,betaFECI))
  colnames(pop$betaFE) = MCMC_Obj$names$FE
  rownames(pop$betaFE)[1] = 'mean'



  # post processing of the profile clustering and outcome effect
  if(comp_cooc){
    cooc = create_co_occurrence_matrix_cpp(MCMC_Obj$Z)
    if (modeClus=='LS'){
      idx = find_ls_optimal_partition(cooc,
                                      MCMC_Obj$Z)
      Zstar = MCMC_Obj$Z[,idx]
      Kstar = length(unique(Zstar))
      mode = 'LS'
      clus = T
    }else if (modeClus=='NG'){
      tmp = estimate_k(cooc, maxk = 15,showplots=F)
      diffs <- diff(tmp$evals)
      diffs <- diffs[-1]
      Kstar <- which.max(abs(diffs[1:15 - 1])) + 1
      Zstar = cluster_similarity(cooc, k = Kstar, specalg = "Ng")

      mode = 'NG'
      clus = T
    }else{
      print('No valid clustering method provided returning only the coocurency
          matrix')
      rep_clust = NULL
      clus = F
    }

    if(clus){
      # Compute the cluster caracteristics
      pvecPost = array(0,dim =c(nUCat,Kstar,nSim))
      centroids = array(0,dim =c(nUCont,Kstar,nSim))
      coVar = array(0,dim = c(nUCont,nUCont,Kstar))
      gamma = array(0,dim =c(qLat,Kstar,nSim))
      if (!is.null(Zstar)){
        optTransf = Zstar*0
        j = 1
        for (it in 1:nSim){
          cIdx = 1
          for (c_ in unique(Zstar)){
            pvecPost[,cIdx,j] = apply( MCMC_Obj$pvec[, MCMC_Obj$Z[,it][which(Zstar == c_)]+1,it,drop=F],1,mean)
            centroids[,cIdx,j] = apply( MCMC_Obj$muClus[, MCMC_Obj$Z[,it][which(Zstar == c_)]+1,it,drop=F],1,mean)
            coVar[,,cIdx] =  coVar[,,cIdx] +apply( MCMC_Obj$PhiClus[[it]][,,MCMC_Obj$Z[,it][which(Zstar == c_)]+1,drop=F],c(1, 2), mean)/nSim
            gamma[,cIdx,j] = apply(MCMC_Obj$gamma[,MCMC_Obj$Z[,it][which(Zstar == c_)]+1,it,drop=F],1,mean)
            optTransf[Zstar == c_] =  cIdx
            cIdx = cIdx+1
          }
          j=j+1
        }

        cen = array(apply(centroids, 2, rowMeans), dim = c(dim(centroids)[1], dim(centroids)[2]))
        pvec = array(apply(pvecPost, 2, rowMeans), dim = c(dim(pvecPost)[1], dim(pvecPost)[2]))
        gam = array(apply(gamma, 2, rowMeans), dim = c(dim(gamma)[1], dim(gamma)[2]))
      }
      rep_clust = list(Zstar = Zstar,
                       Kstar = Kstar,
                       mode = mode,
                       cenStar = centroids,
                       pvecStar = pvecPost,
                       gammaStar = gamma,
                       coVar = coVar,
                       cen = cen,
                       gamma = gam,
                       pvec = pvec)}
  }else{
    rep_clust = NULL
  }
  return(list(
    coocMat = cooc,
    clust = rep_clust,
    pop = pop
  ))
}




