#include "dataClass.h"
// DataObj implementation

DataObj::DataObj(arma::vec Y_, arma::mat XFE_, arma::mat XRE_, arma::mat XL_, arma::mat U_, arma::vec ZRE_, int qRE_ , int nC_) {
  Y = Y_;
  XRE = XRE_;
  XFE = XFE_;
  U = U_;
  XL = XL_;
  ZRE = ZRE_;
  n = Y.n_elem;
  qFE = XFE_.n_cols;
  qRE = qRE_;
  qL = XL_.n_cols;
  qU = U_.n_cols;
  if(qRE>0){
  nRE = ZRE_.max() + 1;}else{
    nRE = 0;
  }
  nC = nC_;
}
