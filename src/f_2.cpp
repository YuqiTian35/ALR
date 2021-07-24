// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
/* # cpp function for elements in the second estimating equation */
List f_2(NumericVector y, int n, NumericMatrix z, NumericVector mu, NumericVector gamma, NumericVector theta, int ncol_z){
  int len = mu.size();
  int n_choose_2 = n * (n - 1) / 2;
  float zeta;
  NumericMatrix inv_s(n_choose_2);
  NumericMatrix s(n_choose_2);
  NumericVector r(n_choose_2);
  NumericMatrix t(n_choose_2, ncol_z);
  float f;
  float nu;
  float offset;
  float temp;
  List res;
  
  int l = 0;
  for (int i = 0; i < len-1; ++i){
    for (int j = i+1; j < len; ++j){
      f = 1 - (1 - theta(l)) * (mu(i) + mu(j));
      nu = (f - std::sqrt(std::pow(f, 2) - 4 * theta(l) * (theta(l) - 1) * mu(i) * mu(j))) / (2 * (theta(l) - 1));
      offset = std::log((mu(i) - nu) / (1 - mu(i) - mu(j) + nu));
      temp = arma::as_scalar(gamma(l) * y(j)) + offset;
      zeta = 1 / (1 + std::exp(-temp));
      inv_s(l, l) = 1 / (zeta * (1 - zeta));
      r(l) = y(i) - zeta;
      t(l, _) = (std::exp(-temp) / (std::pow(1 + std::exp(-temp), 2))) * z(l, _) * y(j);
      ++l;
    }
  }
  
  res["inv_s"] = inv_s;
  res["r"] = r;
  res["t"] = t;
  return res;
}




