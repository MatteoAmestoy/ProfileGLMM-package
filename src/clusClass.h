#pragma once

#include <RcppArmadillo.h>
#include "UtilityFunctions.h" // For helper functions
#include "dataClass.h"

// Parameters for the clustering component
class ParamClus {
public:
  double lam0;
  arma::vec mu0;
  double nu0;
  arma::mat Psi0;
  arma::mat mu;
  arma::cube Sigma;
  ParamClus(double lam0, arma::vec mu0, double nu0, arma::mat Psi0, arma::mat mu, arma::cube Sigma);
  void update(DataObj data, arma::ivec Z);
};
