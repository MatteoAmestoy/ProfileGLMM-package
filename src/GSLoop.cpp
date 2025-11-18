#include <RcppArmadillo.h>
#include <mvnorm.h>
#include <wishart.h>

#include "assignClass.h"
#include "clusClass.h"
#include "dataClass.h"
#include "GLMMClass.h"
#include "UtilityFunctions.h"

using namespace Rcpp;

// [[Rcpp::export]]
List GSLoopCPP(
    int nIt,
    int nBurnIn,
    int nC,
    int qRE,
    int qUCont, //---------------------
    arma::vec Y,
    arma::mat XFE,
    arma::mat XRE,
    arma::mat XL,
    arma::mat UCont,//---------------------
    arma::mat UCat,//---------------------
    arma::vec catInd,//---------------------
    arma::vec ZRE,
    arma::vec beta,
    double sig2,
    arma::mat WRE,
    arma::mat muClus,
    arma::cube SigmaClus,
    arma::mat pvecClus,//---------------------
    arma::mat WLat,
    arma::mat gammaLat,
    double a,
    double b,
    double lambdaFE,
    arma::mat PhiRE,
    double etaRE,
    double lam0,
    arma::vec mu0,
    double nu0,
    arma::mat Psi0,
    arma::vec alpha0,//---------------------
    arma::mat PhiLat,
    double etaLat,
    double scale,
    double shape,
    int regType) {
  // Data object
  DataObj data(Y, XFE, XRE, XL, UCont, UCat, catInd, ZRE, qRE, qUCont, nC);

  // Initialize theta
  ParamGLMM thetaLMM(beta, sig2, gammaLat, WLat, WRE, a, b, lambdaFE, PhiLat,
                     etaLat, PhiRE, etaRE, data.XFE);
  ParamClus thetaClus(lam0, mu0, nu0, Psi0, muClus, SigmaClus,
                      alpha0, pvecClus);

  arma::vec p0(nC);
  arma::ivec Z(data.n);
  Z = arma::randi(data.n, arma::distr_param(0, nC - 1));

  ParamAssign thetaAssign(Z, p0, scale, shape);

  int nStore = nIt - nBurnIn - 1;
  arma::cube StoreWL(data.qL, data.qL, nStore);

  arma::cube StoreWRE(std::max(1,data.qRE), std::max(1,data.qRE), nStore);

  arma::imat StoreZ(data.n, nStore);
  arma::cube StoremuClus(data.qUCont, nC, nStore);
  arma::mat StoreBeta(data.qFE, nStore);
  arma::vec Storesig2(nStore);
  arma::vec Storealpha(nStore);
  arma::cube StoreGamma(data.qL, nC, nStore);
  arma::cube Storepvec(thetaClus.pVec.n_rows, thetaClus.pVec.n_cols, nStore);
  arma::field<arma::cube> StorePhiClus(nStore);
  int storeIdx = 0;

  if (regType == 0) {
    for (int it = 0; it < nIt; it++) {
      thetaLMM.updateLinear(data, thetaAssign.Z, thetaAssign.cluster_count);
      thetaClus.update(data, thetaAssign.Z);
      thetaAssign.updateLinear(data, data.Y - thetaLMM.YRE - thetaLMM.YFE, thetaLMM.sig2,
                               thetaClus.Sigma,
                               thetaClus.mu,
                               thetaLMM.gamma,
                               thetaClus.pVec);
      if ((it % 1000) == 0) {
        Rcout << "Iteration: " << it << "\n";
        Rcout << "Nb of non 0 clusters: " << thetaAssign.non_0_clust << "\n";
      }
      if (it > nBurnIn) {
        StoreWL.slice(storeIdx) = thetaLMM.WLat;
        if(data.qRE>0){
          StoreWRE.slice(storeIdx) = thetaLMM.WRE;}
        StorePhiClus(storeIdx) = thetaClus.Sigma;
        StoreZ.col(storeIdx) = thetaAssign.Z;
        StoremuClus.slice(storeIdx) = thetaClus.mu;
        StoreBeta.col(storeIdx) = thetaLMM.beta;
        StoreGamma.slice(storeIdx) = thetaLMM.gamma;
        Storesig2(storeIdx) = thetaLMM.sig2;
        Storealpha(storeIdx) = thetaAssign.alpha;
        Storepvec.slice(storeIdx) = thetaClus.pVec;
        storeIdx += 1;
      }
    }
  } else if (regType == 1) {
    for (int it = 0; it < nIt; it++) {
      thetaLMM.updateProbit(data, thetaAssign.Z, thetaAssign.cluster_count);
      thetaClus.update(data, thetaAssign.Z);
      thetaAssign.updateProbit(data, thetaLMM.YRE, thetaLMM.YFE,
                               thetaClus.Sigma,
                               thetaClus.mu,
                               thetaLMM.gamma,
                               thetaClus.pVec);
      if ((it % 1000) == 0) {
        Rcout << "Iteration: " << it << "\n";
        Rcout << "Nb of non 0 clusters: " << thetaAssign.non_0_clust << "\n";
      }
      if (it > nBurnIn) {
        StoreWL.slice(storeIdx) = thetaLMM.WLat;
        if(data.qRE>0){
          StoreWRE.slice(storeIdx) = thetaLMM.WRE;}
        StorePhiClus(storeIdx) = thetaClus.Sigma;
        StoreZ.col(storeIdx) = thetaAssign.Z;
        StoremuClus.slice(storeIdx) = thetaClus.mu;
        StoreBeta.col(storeIdx) = thetaLMM.beta;
        StoreGamma.slice(storeIdx) = thetaLMM.gamma;
        Storesig2(storeIdx) = thetaLMM.sig2;
        Storealpha(storeIdx) = thetaAssign.alpha;
        Storepvec.slice(storeIdx) = thetaClus.pVec;
        storeIdx += 1;
      }
    }
  }

  List Store = List::create(Named("WL") = StoreWL,
                            _["WRE"] = StoreWRE,
                            _["Z"] = StoreZ,
                            _["muClus"] = StoremuClus,
                            _["PhiClus"] = StorePhiClus,
                            _["beta"] = StoreBeta,
                            _["sig2"] = Storesig2,
                            _["zeta"] = Storealpha,
                            _["pvec"] = Storepvec,
                            _["gamma"] = StoreGamma);

  return Store;
}
