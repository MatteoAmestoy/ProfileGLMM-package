#pragma once

#include <RcppArmadillo.h>
#include "UtilityFunctions.h" // For helper functions
#include "dataClass.h"

// Parameters for the assignment component
class ParamAssign {
public:
  arma::ivec Z;
  arma::vec p0;
  double alpha;
  double scale;
  double shape;
  arma::vec cluster_count;
  int non_0_clust;
  arma::sp_umat adjMat;
  ParamAssign(arma::ivec Z, arma::vec p0, double scale, double shape);
  void updateLinear(DataObj data, arma::vec Y, double sig2, arma::cube SigmaGM,
                    arma::mat muGM, arma::mat gammaL, arma::mat pvec);
  void updateProbit(DataObj data, arma::vec YFE, arma::vec YRE, arma::cube SigmaGM,
                    arma::mat muGM, arma::mat gammaL, arma::mat pvec);
};
