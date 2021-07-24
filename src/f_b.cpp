// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
/* # cpp function for b in `func_b */
NumericVector f_b(NumericVector nu, NumericVector mu){
  int len = mu.size();
  
  NumericMatrix b(len);
  for (int i=0; i<len; i++){
    b(i,i) = mu(i) * (1 - mu(i));
  }
  
  int l = 0;
  for (int i = 0; i < len-1; ++i){
    for (int j = i+1; j < len; ++j){
      b(i,j) = nu(l) - mu(i) * mu(j);
      b(j,i) = b(i,j);
      ++l;
    }
  }
  return b;
}