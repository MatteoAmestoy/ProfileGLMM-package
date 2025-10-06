#include "GLMMClass.h"
#include <wishart.h> // Assuming this is an external header for wishart distribution
#include <truncnorm.h>

// ParamGLMM implementation
ParamGLMM::ParamGLMM(arma::vec beta_, double sig2_, arma::mat gamma_, arma::mat WLat_, arma::mat WRE_, double a_, double b_, double lambda_, arma::mat PhiLat_, double etaLat_, arma::mat PhiRE_, double etaRE_, arma::mat XFE) {
  int n = XFE.n_rows;
  beta = beta_;
  gamma = gamma_;
  WLat = WLat_;
  WRE = WRE_;
  a = a_;
  b = b_;
  lambda = lambda_;
  PhiLat = PhiLat_;
  etaLat = etaLat_;
  PhiRE = PhiRE_;
  etaRE = etaRE_;
  arma::vec YLat_(n);
  arma::vec YRE_(n);
  arma::vec prob_score_(n);
  YLat = YLat_;
  YRE = YRE_;
  prob_score = prob_score_;
  YFE = YRE_;
  sig2 = sig2_;
}

void ParamGLMM::updateLinear(DataObj data, arma::ivec Z, arma::vec cluster_count) {
  // We start by jointly updating the latent and fixed effect parts of the model and then do the RE
  arma::vec Y = data.Y - YRE;
  arma::mat W_m = arma::inv(WLat);
  arma::mat omega = PhiLat;
  arma::mat B = arma::eye(data.qFE, data.qFE);
  arma::vec bhat(data.qFE);
  B = B * lambda / sig2;
  arma::cube Plat(data.qL, data.qFE, data.nC);
  for(int c = 0; c < data.nC; c++) {
    if(cluster_count(c) > 0) {
      arma::uvec idx = arma::find(Z == c);
      arma::mat XL_ = data.XL.rows(idx);
      arma::mat X_ = data.XFE.rows(idx);
      arma::mat XXL = X_.t() * XL_;
      arma::mat XLXL = XL_.t() * XL_;
      arma::mat VmCore = arma::inv(W_m + XLXL / sig2) / pow(sig2, 2);
      arma::mat XLY = XL_.t() * Y(idx);
      B = B + (X_.t() * X_ / sig2 - XXL * VmCore * XXL.t());
      bhat = bhat + (X_.t() * Y(idx) / sig2 - XXL * VmCore * XLY);
      arma::mat Sig = WLat - WLat * XLXL * WLat / sig2 + WLat * XLXL * VmCore * XLXL * WLat;
      gamma.col(c) = mvrnormArma(1, WLat * XLY / sig2 - WLat * XLXL * VmCore * XLY, (Sig + Sig.t()) / 2.0).t();
      Plat.slice(c) = WLat * XXL.t() / sig2 - WLat * XLXL * VmCore * XXL.t();
    }
  }
  B = arma::inv(B);
  bhat = B * bhat;
  beta = mvrnormArma(1, bhat, (B + B.t()) / 2.0).t();
  YFE = data.XFE * beta;
  for(int c = 0; c < data.nC; c++) {
    if(cluster_count(c) > 0) {
      gamma.col(c) = gamma.col(c) - Plat.slice(c) * beta;
      arma::uvec idx = arma::find(Z == c);
      YLat(idx) = data.XL.rows(idx) * gamma.col(c);
    } else {
      arma::vec mu(data.qL);
      gamma.col(c) = mvrnormArma(1, mu, (WLat + WLat.t()) / 2.0).t();
    }
    omega += gamma.col(c) * gamma.col(c).t();
  }
  WLat = riwish(data.nC + etaLat, (omega + omega.t()) / 2.0);
  double b_ = b + arma::as_scalar((Y - YFE - YLat).t() * (Y - YFE - YLat)) / 2.0;
  sig2 = 1.0 / (Rcpp::rgamma(1, a + data.n / 2.0, 1 / b_)(0));

  // Move to Random effects
  if(data.qRE>0){
    Y = data.Y - YFE - YLat;
    arma::mat gammaRE(data.qRE, data.nRE);
    arma::mat omegaRE = PhiRE;
    arma::mat WRE_m = inv(WRE);
    for(int ind = 0; ind < data.nRE; ind++) {
      arma::uvec idxRE = arma::find(data.ZRE == ind);
      arma::mat XRE_ = data.XRE.rows(idxRE);
      arma::mat XX = XRE_.t() * XRE_;
      arma::mat VmCore2RE = WRE * XX * arma::inv(WRE_m + XX / sig2) / pow(sig2, 2);
      arma::mat XY = XRE_.t() * Y(idxRE);
      arma::mat SigRE = WRE - WRE * XX * WRE / sig2 + VmCore2RE * XX * WRE;
      gammaRE.col(ind) = mvrnormArma(1, WRE * XY / sig2 - VmCore2RE * XY, (SigRE + SigRE.t()) / 2).t();
      omegaRE += gammaRE.col(ind) * gammaRE.col(ind).t();
      YRE(idxRE) = XRE_ * gammaRE.col(ind);
    }
    WRE = riwish(data.nRE + etaRE, omegaRE);}
}

void ParamGLMM::updateProbit(DataObj data, arma::ivec Z, arma::vec cluster_count) {
  // Probit update logic
  double amin;
  double amax;
  arma::vec mu_ = YRE+YFE+YLat;
  for(int o = 0; o < data.n; o++){
    amin = -pow(10,6)*(1-data.Y(o));
    amax = pow(10,6)*data.Y(o);
    prob_score(o) = r_truncnorm(mu_(o), 1, amin,
               amax);
  }
  arma::vec Y = prob_score - YRE;
  arma::mat W_m = arma::inv(WLat);
  arma::mat omega = PhiLat;
  arma::mat B = arma::eye(data.qFE, data.qFE);
  arma::vec bhat(data.qFE);
  B = B * lambda;
  arma::cube Plat(data.qL, data.qFE, data.nC);
  for(int c = 0; c < data.nC; c++) {
    if(cluster_count(c) > 0) {
      arma::uvec idx = arma::find(Z == c);
      arma::mat XL_ = data.XL.rows(idx);
      arma::mat X_ = data.XFE.rows(idx);
      arma::mat XXL = X_.t() * XL_;
      arma::mat XLXL = XL_.t() * XL_;
      arma::mat VmCore = arma::inv(W_m + XLXL) ;
      arma::mat XLY = XL_.t() * Y(idx);
      B = B + (X_.t() * X_  - XXL * VmCore * XXL.t());
      bhat = bhat + (X_.t() * Y(idx)  - XXL * VmCore * XLY);
      arma::mat Sig = WLat - WLat * XLXL * WLat  + WLat * XLXL * VmCore * XLXL * WLat;
      gamma.col(c) = mvrnormArma(1, WLat * XLY  - WLat * XLXL * VmCore * XLY, (Sig + Sig.t()) / 2.0).t();
      Plat.slice(c) = WLat * XXL.t() - WLat * XLXL * VmCore * XXL.t();
    }
  }
  B = arma::inv(B);
  bhat = B * bhat;
  beta = mvrnormArma(1, bhat, (B + B.t()) / 2.0).t();
  YFE = data.XFE * beta;
  for(int c = 0; c < data.nC; c++) {
    if(cluster_count(c) > 0) {
      gamma.col(c) = gamma.col(c) - Plat.slice(c) * beta;
      arma::uvec idx = arma::find(Z == c);
      YLat(idx) = data.XL.rows(idx) * gamma.col(c);
    } else {
      arma::vec mu(data.qL);
      gamma.col(c) = mvrnormArma(1, mu, (WLat + WLat.t()) / 2.0).t();
    }
    omega += gamma.col(c) * gamma.col(c).t();
  }
  WLat = riwish(data.nC + etaLat, (omega + omega.t()) / 2.0);

  // mu_ = YRE+YFE+YLat;
  // for(int o = 0; o < data.n; o++){
  //   amin = -pow(10,6)*(1-data.Y(o));
  //   amax = pow(10,6)*data.Y(o);
  //   prob_score(o) = r_truncnorm(mu_(o), 1, amin,
  //              amax);
  // }


  if(data.qRE>0){
    Y = prob_score - YFE - YLat;
    arma::mat gammaRE(data.qRE, data.nRE);
    arma::mat omegaRE = PhiRE;
    arma::mat WRE_m = inv(WRE);
    for(int ind = 0; ind < data.nRE; ind++) {
      arma::uvec idxRE = arma::find(data.ZRE == ind);
      arma::mat XRE_ = data.XRE.rows(idxRE);
      arma::mat XX = XRE_.t() * XRE_;
      arma::mat VmCore2RE = WRE * XX * arma::inv(WRE_m + XX );
      arma::mat XY = XRE_.t() * Y(idxRE);
      arma::mat SigRE = WRE - WRE * XX * WRE+ VmCore2RE * XX * WRE;
      gammaRE.col(ind) = mvrnormArma(1, WRE * XY - VmCore2RE * XY, (SigRE + SigRE.t()) / 2).t();
      omegaRE += gammaRE.col(ind) * gammaRE.col(ind).t();
      YRE(idxRE) = XRE_ * gammaRE.col(ind);
    }
    WRE = riwish(data.nRE + etaRE, omegaRE);}

}
