#pragma once

#include <RcppArmadillo.h>
#include "UtilityFunctions.h" // For helper functions
#include "dataClass.h"

using namespace Rcpp;

// Parameter container for the GLMM
class ParamGLMM {
public:
  arma::vec beta;
  double sig2;
  arma::mat gamma;
  arma::mat WLat;
  arma::mat WRE;
  double a;
  double b;
  double lambda;
  arma::mat PhiLat;
  double etaLat;
  arma::mat PhiRE;
  double etaRE;
  arma::vec prob_score;
  arma::vec YFE;
  arma::vec YRE;
  arma::vec YLat;
  ParamGLMM(arma::vec beta, double sig2, arma::mat gamma, arma::mat WLat, arma::mat WRE, double a, double b, double lambda, arma::mat PhiLat, double etaLat, arma::mat PhiRE, double etaRE, arma::mat XFE);
  void updateLinear(DataObj data, arma::ivec Z, arma::vec cluster_count);
  void updateProbit(DataObj data, arma::ivec Z, arma::vec cluster_count);
};
