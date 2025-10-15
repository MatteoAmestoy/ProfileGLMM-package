#include "clusClass.h"
#include <wishart.h> // Assuming this is an external header for wishart distribution

// ParamClus implementation
ParamClus::ParamClus(double lam0_, arma::vec mu0_, double nu0_, arma::mat Psi0_, arma::mat mu_, arma::cube Sigma_,
                     arma::vec alpha0_,arma::mat pVec_) {
  lam0 = lam0_;
  nu0 = nu0_;
  mu0 = mu0_;
  Psi0 = Psi0_;
  mu = mu_;
  Sigma = Sigma_;
  pVec = pVec_;
  alpha0 = alpha0_;
}

void ParamClus::update(DataObj data, arma::ivec Z) {
  arma::mat EUc(data.qUCont, data.nC);
  arma::vec nCvec(data.nC);
  arma::cube vUc(data.qUCont, data.qUCont, data.nC);
  EUc.fill(0);
  vUc.fill(0);


  if (data.UContBool){
    for(int i = 0; i < data.n; i++) {
      EUc.col(Z(i)) += data.UCont.row(i).t();
      nCvec(Z(i)) += 1.0;
      vUc.slice(Z(i)) += data.UCont.row(i).t() * data.UCont.row(i);
    }

    for(int c = 0; c < data.nC; c++) {
      arma::vec mun = mu0;
      double lambdan = lam0;
      double nun = nu0;
      arma::mat Psin = Psi0;

      if (nCvec(c) != 0) {
        mun = (lam0 * mu0 + EUc.col(c)) / (lam0 + nCvec(c));
        lambdan += nCvec(c);
        nun += nCvec(c);
        Psin += ((EUc.col(c) / nCvec(c) - mu0) * (EUc.col(c) / nCvec(c) - mu0).t()) * (lam0 * nCvec(c)) / (lam0 + nCvec[c]) +
          vUc.slice(c) - EUc.col(c) * EUc.col(c).t() / nCvec(c);
      }
      Sigma.slice(c) = riwish(nun, (Psin + Psin.t()) / 2.0);
      mu.col(c) = mvrnormArma(1, mun, (Sigma.slice(c) + Sigma.slice(c).t()) / 2.0 / lambdan).t();
    }
  }
  if (data.UCatBool){

    //Build the current count per cluster
    arma::mat p(data.UCat.n_cols,data.nC);
    for(int c = 0; c < data.nC; c++) {
      arma::uvec  idx = arma::find(Z == c);
      p.col(c) = sum(data.UCat.rows(idx),0);
    }


    for(int c = 0; c < data.nC; c++) {
      arma::vec p_ = p.col(c);
      for (int cat = 0; cat < data.nCat; cat++) {
        arma::uvec idxCat = arma::find(data.catInd == cat);
        p_(idxCat) = rdirichlet_cpp(p_(idxCat)+alpha0(cat)*arma::ones(data.nUCat(cat)));
      }
      pVec.col(c) = p_;
    }
  }
}
