#pragma once

#include <RcppArmadillo.h>
#include "UtilityFunctions.h" // For helper functions

using namespace Rcpp;

// Data container for the GLMM
class DataObj {
public:
  arma::vec Y; // outcome
  arma::mat XFE;// Fixed Effects covariates
  arma::mat XRE;// Random Effects covariates
  arma::mat XL;// Interaction Effects covariates
  arma::mat UCont;// Continuous clustering covariates
  ListMatrix UCat;// Categorical clustering covariates
  arma::vec ZRE; // Random effect statisical unit membership
  arma::mat PX;
  arma::mat VmX;
  int n;// nb observations
  int qFE;// nb FE covariates
  int qRE;// nb RE covariates
  int qL;// nb Interaction covariates
  int nRE;// nb RE stat units
  int qUCont;//nb continuous clustering covariates
  int nC;//max nb of latent clusters
  int nUCat; //nb of categorical covariates
  bool UCatBool;//presence or absence of continuous clustering covariates
  bool UContBool;//presence or absence of cat clustering covariates

  DataObj(arma::vec Y, arma::mat XFE, arma::mat XRE, arma::mat XL, arma::mat U, arma::vec ZRE, int qRE, int nC);
};
