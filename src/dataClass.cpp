#include "dataClass.h"
// DataObj implementation

DataObj::DataObj(arma::vec Y_, arma::mat XFE_, arma::mat XRE_, arma::mat XL_, arma::mat UCont_,
                 arma::mat UCat_, arma::vec catInd_, arma::vec ZRE_, int qRE_ ,
                 int qUCont_ , int nC_) {
  Y = Y_;
  XRE = XRE_;
  XFE = XFE_;
  catInd = catInd_;
  nCat = catInd_.max()+1;
  UCatBool  = false;
  if (nCat>0){
    UCatBool= true;
    UCat = UCat_;
    for( int cat = 0; cat < nCat; cat++){
      arma::uvec idx = arma::find(catInd_ == cat);
      nUCat(cat) = idx.n_elem;
    }
  }
  XL = XL_;
  ZRE = ZRE_;
  n = Y.n_elem;
  qFE = XFE_.n_cols;
  qRE = qRE_;
  qL = XL_.n_cols;
  UContBool = false;
  qUCont = qUCont_;
  if (qUCont>0){
    UCont = UCont_;
    UContBool = true;
  }

  if(qRE>0){
    nRE = ZRE_.max() + 1;}else{
      nRE = 0;
    }
    nC = nC_;
}
