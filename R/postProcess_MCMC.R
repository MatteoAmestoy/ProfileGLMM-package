#' Post processing of the MCMC chain.
#'
#' @param MCMCObj Profile GLMM  MCMC output of  profileGLMM_Gibbs function
#' @returns The Gibbs-sampled posterior. Output of the cpp GSLoopCPP function.
#' @export
#'
#' @examples
profileGLMM_postProcess = function(MCMCObj,modeClus='NG'){

  nSim = dim(MCMCObj$Z)[2]
  nUCat = dim(MCMCObj$pvec)[1]
  nUCont = dim(MCMCObj$muClus)[1]
  qLat = dim(MCMCObj$gamma)[1]


  cooc = create_co_occurrence_matrix_cpp(MCMCObj$Z)
  if (modeClus=='LS'){
    idx = find_ls_optimal_partition(cooc,
                                    MCMCObj$Z)
    Zstar = MCMCObj$Z[,idx]
    Kstar = length(unique(Zstar))
    mode = 'LS'
  }else if (modeClus=='NG'){
    tmp = estimate_k(cooc, maxk = 15,showplots=F)
    diffs <- diff(tmp$evals)
    diffs <- diffs[-1]
    Kstar <- which.max(abs(diffs[1:15 - 1])) + 1
    Zstar = cluster_similarity(cooc, k = Kstar, specalg = "Ng")

    mode = 'NG'
  }

  # Compute the cluster caracteristics
  pvecPost = array(0,dim =c(nUCat,Kstar,nSim))
  centroids = array(0,dim =c(nUCont,Kstar,nSim))
  coVar = array(0,dim = c(nUCont,nUCont,Kstar))
  gamma = array(0,dim =c(qLat,Kstar,nSim))

  optTransf = Zstar*0
  j = 1
  for (it in 1:nSim){
    cIdx = 1
    for (c_ in unique(Zstar)){
      pvecPost[,cIdx,j] = apply( MCMCObj$pvec[, MCMCObj$Z[,it][which(Zstar == c_)]+1,it,drop=F],1,mean)
      centroids[,cIdx,j] = apply( MCMCObj$muClus[, MCMCObj$Z[,it][which(Zstar == c_)]+1,it,drop=F],1,median)
      coVar[,,cIdx] =  coVar[,,cIdx] +median( MCMCObj$PhiClus[[it]][,,MCMCObj$Z[,it][which(Zstar == c_)]+1])/nSim
      gamma[,cIdx,j] = apply(MCMCObj$gamma[,MCMCObj$Z[,it][which(Zstar == c_)]+1,it,drop=F],1,mean)
      optTransf[Zstar == c_] =  cIdx
      cIdx = cIdx+1
    }
    j=j+1
  }
  cen = apply(centroids, 2, rowMeans)
  pvec = apply(pvecPost, 2, rowMeans)
  gam = apply(gamma, 2, rowMeans)

  return(list(
    coocMat = cooc,
    Zstar = Zstar,
    Kstar = Kstar,
    mode = mode,
    cenStar = centroids,
    pvecStar = pvecPost,
    gammaStar = gamma,
    coVar = coVar,
    cen = cen,
    gamma = gam,
    pvec = pvec

  ))
}
