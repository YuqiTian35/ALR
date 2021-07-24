// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat mat_inv(arma::mat a) { return inv(a); }
 

// [[Rcpp::export]] 
arma::mat mat_3(arma::mat a, arma::mat b, arma::mat c) { return a * b * c; }
    

// [[Rcpp::export]]
arma::mat mat_sqr(arma::mat a) { return a * a.t(); }
