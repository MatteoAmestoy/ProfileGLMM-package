#include "assignClass.h"


// ParamAssign implementation
ParamAssign::ParamAssign(arma::ivec Z_, arma::vec p0_, double scale_, double shape_) {
  Z = Z_;
  p0 = p0_;
  int nC = p0.n_elem;
  int n = Z.n_elem;
  alpha = R::rgamma(shape_, scale_);
  scale = scale_;
  shape = shape_;
  arma::vec cluster_count_(nC);
  for(int i = 0; i < n; i++) {
    cluster_count_(Z(i)) += 1;
  }
  cluster_count = cluster_count_;
  int non_0_clust_loc = 0;
  for(int c = 0; c < nC; c++) {
    if (cluster_count(c) > 0) {
      non_0_clust_loc += 1.0;
    }
  }
  non_0_clust = non_0_clust_loc;
  adjMat = arma::sp_umat(n, n);
}

void ParamAssign::updateLinear(DataObj data, arma::vec Y, double sig2, arma::cube SigmaGM, arma::mat muGM, arma::mat gammaL) {
  // Update logic for linear model
  double VarAcc;
  arma::vec V(data.nC);
  VarAcc = 0;
  for(int c = 0; c < data.nC; c++) {
    V(data.nC - 1 - c) = R::rbeta(1 + cluster_count(data.nC - 1 - c), alpha + VarAcc);
    VarAcc += cluster_count(data.nC - 1 - c);
  }
  VarAcc = 1;
  for(int c = 0; c < data.nC; c++) {
    p0(c) = V(c) * VarAcc;
    VarAcc = VarAcc * (1 - V(c));
  }
  arma::vec lp0(data.nC);
  arma::cube SigmaGMinvC(data.qU, data.qU, data.nC);
  arma::mat predC(data.n, data.nC);
  arma::vec cluster_count_loc(data.nC);

  for(int c = 0; c < data.nC; c++) {
    lp0(c) = log(p0(c)) - arma::log_det_sympd(SigmaGM.slice(c)) / 2;
    SigmaGMinvC.slice(c) = arma::inv(SigmaGM.slice(c));
    predC.col(c) = Y - data.XL * gammaL.col(c);
  }
  arma::vec logprob(data.nC);
  for(int o = 0; o < data.n; o++) {
    logprob = lp0;
    double nconst = 0;
    for(int c = 0; c < data.nC; c++) {
      logprob(c) += -arma::as_scalar((data.U.row(o) - muGM.col(c).t()) * SigmaGMinvC.slice(c) * (data.U.row(o).t() - muGM.col(c))) / 2.0;
      logprob(c) += -pow(predC(o, c), 2) / sig2 / 2.0;
      logprob(c) = exp(logprob(c));
      nconst += logprob(c);
    }
    double u = runif(1, 0, 1)(0) * nconst;
    for(int c = 0; c < data.nC; c++) {
      if (u < logprob(c)) {
        Z(o) = c;
        cluster_count_loc(c) += 1.0;
        u = nconst;
      } else {
        u = u - logprob(c);
      }
    }
  }
  int non_0_clust_loc = 0;
  for(int c = 0; c < data.nC; c++) {
    if (cluster_count_loc(c) > 0) {
      non_0_clust_loc += 1.0;
    }
  }
  alpha = R::rgamma(shape + data.nC, 1 / (1 / scale - log(VarAcc)));
  cluster_count = cluster_count_loc;
  non_0_clust = non_0_clust_loc;
}

void ParamAssign::updateProbit(DataObj data, arma::vec YFE, arma::vec YRE, arma::cube SigmaGM, arma::mat muGM, arma::mat gammaL) {
  // Probit update logic
  // Update logic for linear model
  double VarAcc;
  arma::vec V(data.nC);
  VarAcc = 0;
  for(int c = 0; c < data.nC; c++) {
    V(data.nC - 1 - c) = R::rbeta(1 + cluster_count(data.nC - 1 - c), alpha + VarAcc);
    VarAcc += cluster_count(data.nC - 1 - c);
  }
  VarAcc = 1;
  for(int c = 0; c < data.nC; c++) {
    p0(c) = V(c) * VarAcc;
    VarAcc = VarAcc * (1 - V(c));
  }
  arma::vec lp0(data.nC);
  arma::cube SigmaGMinvC(data.qU, data.qU, data.nC);
  arma::mat predC(data.n, data.nC); // predicted mean for each cluster based on cluster membership
  arma::vec cluster_count_loc(data.nC);

  for(int c = 0; c < data.nC; c++) {
    lp0(c) = log(p0(c)) - arma::log_det_sympd(SigmaGM.slice(c)) / 2;
    SigmaGMinvC.slice(c) = arma::inv(SigmaGM.slice(c));
    predC.col(c) = arma::normcdf(YFE + YRE + data.XL * gammaL.col(c),0,1);//prob that Y=1
    predC.col(c) = predC.col(c)%data.Y+(1-data.Y)%(1-predC.col(c));
  }
  arma::vec logprob(data.nC);
  for(int o = 0; o < data.n; o++) {
    logprob = lp0;
    double nconst = 0;
    for(int c = 0; c < data.nC; c++) {
      logprob(c) += -arma::as_scalar((data.U.row(o) - muGM.col(c).t()) * SigmaGMinvC.slice(c) * (data.U.row(o).t() - muGM.col(c))) / 2.0;
      logprob(c) = exp(logprob(c));
      logprob(c) = logprob(c)*predC(o,c);
      nconst += logprob(c);
    }
    double u = runif(1, 0, 1)(0) * nconst;
    for(int c = 0; c < data.nC; c++) {
      if (u < logprob(c)) {
        Z(o) = c;
        cluster_count_loc(c) += 1.0;
        u = nconst;
      } else {
        u = u - logprob(c);
      }
    }
  }
  int non_0_clust_loc = 0;
  for(int c = 0; c < data.nC; c++) {
    if (cluster_count_loc(c) > 0) {
      non_0_clust_loc += 1.0;
    }
  }
  alpha = R::rgamma(shape + data.nC, 1 / (1 / scale - log(VarAcc)));
  cluster_count = cluster_count_loc;
  non_0_clust = non_0_clust_loc;
  }
