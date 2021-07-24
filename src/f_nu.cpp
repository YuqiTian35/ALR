// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
/* # cpp function for nu in `func_nu` */
NumericVector f_nu(NumericVector theta, NumericVector mu){
  int len = mu.size();
  int retLen = len * (len - 1) / 2;
  NumericVector nu(retLen);
  float f;
  
  int l = 0;
  for (int i = 0; i < len -1; ++i){
    for (int j = i+1; j < len; ++j){
      f = 1 - (1 - theta(l)) * (mu(i) + mu(j));
      nu(l) = (f - std::sqrt(std::pow(f, 2) - 4 * theta(l) * (theta(l) - 1) * mu(i) * mu(j))) / (2 * (theta(l) - 1));
      ++l;
    }
  }
  return nu;
};