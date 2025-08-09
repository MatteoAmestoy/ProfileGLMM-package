#include "dataClass.h"
// DataObj implementation

DataObj::DataObj(arma::vec Y_, arma::mat XFE_, arma::mat XRE_, arma::mat XL_, arma::mat U_, arma::vec ZRE_, int nC_) {
  Y = Y_;
  XRE = XRE_;
  XFE = XFE_;
  U = U_;
  XL = XL_;
  ZRE = ZRE_;
  n = Y.n_elem;
  qFE = XFE_.n_cols;
  qRE = XRE_.n_cols;
  qL = XL_.n_cols;
  qU = U_.n_cols;
  nRE = ZRE_.max() + 1;
  nC = nC_;
}