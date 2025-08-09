#pragma once

#include <RcppArmadillo.h>

// Declarations of helper functions
arma::mat rdirichlet_cpp(arma::vec alpha_m);
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);