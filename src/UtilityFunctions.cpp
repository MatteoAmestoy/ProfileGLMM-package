#include "UtilityFunctions.h"
#include <Rcpp.h>

using namespace Rcpp;

//' Draw from a Dirichlet distribution
//'
//' @param alpha_m A vector of parameters.
//' @return A vector of probabilities drawn from a Dirichlet distribution.
arma::mat rdirichlet_cpp(arma::vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  arma::vec distribution = arma::zeros(distribution_size);
  double sum_term = 0;
  for (int j = 0; j < distribution_size; ++j) {
    double cur = R::rgamma(alpha_m[j], 1.0);
    distribution(j) = cur;
    sum_term += cur;
  }
  return(distribution / sum_term);
}

//' Draw from a Multivariate Normal distribution
//'
//' @param n Number of samples.
//' @param mu Mean vector.
//' @param sigma Covariance matrix.
//' @return A matrix with 'n' rows, each a draw from the distribution.
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  bool success = false;
  int breakChol = 0;
  arma::mat R;
  while(success == false and breakChol < 10) {
    success = arma::chol(R, sigma);
    if(success == false) {
      Rcout << "Chol failed adding identity" << "\n";
      breakChol += 1;
      sigma += arma::eye(sigma.n_rows, sigma.n_rows) * 1e-6;
    }
  }
  return arma::repmat(mu, 1, n).t() + Y * R;
}


