#pragma once

#include <RcppArmadillo.h>
#include "UtilityFunctions.h" // For helper functions

using namespace Rcpp;

// Data container for the GLMM
class DataObj {
public:
  arma::vec Y;
  arma::mat XFE;
  arma::mat XRE;
  arma::mat XL;
  arma::mat U;
  arma::vec ZRE;
  arma::mat PX;
  arma::mat VmX;
  int n;
  int qFE;
  int qRE;
  int qL;
  int nRE;
  int qU;
  int nC;
  DataObj(arma::vec Y, arma::mat XFE, arma::mat XRE, arma::mat XL, arma::mat U, arma::vec ZRE, int nC);
};