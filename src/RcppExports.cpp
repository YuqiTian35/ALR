// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// f_2
/* # cpp function for elements in the second estimating equation */ List f_2(NumericVector y, int n, NumericMatrix z, NumericVector mu, NumericVector gamma, NumericVector theta, int ncol_z);
RcppExport SEXP _ALR_f_2(SEXP ySEXP, SEXP nSEXP, SEXP zSEXP, SEXP muSEXP, SEXP gammaSEXP, SEXP thetaSEXP, SEXP ncol_zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type ncol_z(ncol_zSEXP);
    rcpp_result_gen = Rcpp::wrap(f_2(y, n, z, mu, gamma, theta, ncol_z));
    return rcpp_result_gen;
END_RCPP
}
// f_b
/* # cpp function for b in `func_b */ NumericVector f_b(NumericVector nu, NumericVector mu);
RcppExport SEXP _ALR_f_b(SEXP nuSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(f_b(nu, mu));
    return rcpp_result_gen;
END_RCPP
}
// f_nu
/* # cpp function for nu in `func_nu` */ NumericVector f_nu(NumericVector theta, NumericVector mu);
RcppExport SEXP _ALR_f_nu(SEXP thetaSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(f_nu(theta, mu));
    return rcpp_result_gen;
END_RCPP
}
// mat_inv
arma::mat mat_inv(arma::mat a);
RcppExport SEXP _ALR_mat_inv(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_inv(a));
    return rcpp_result_gen;
END_RCPP
}
// mat_3
arma::mat mat_3(arma::mat a, arma::mat b, arma::mat c);
RcppExport SEXP _ALR_mat_3(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_3(a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// mat_sqr
arma::mat mat_sqr(arma::mat a);
RcppExport SEXP _ALR_mat_sqr(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_sqr(a));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ALR_f_2", (DL_FUNC) &_ALR_f_2, 7},
    {"_ALR_f_b", (DL_FUNC) &_ALR_f_b, 2},
    {"_ALR_f_nu", (DL_FUNC) &_ALR_f_nu, 2},
    {"_ALR_mat_inv", (DL_FUNC) &_ALR_mat_inv, 1},
    {"_ALR_mat_3", (DL_FUNC) &_ALR_mat_3, 3},
    {"_ALR_mat_sqr", (DL_FUNC) &_ALR_mat_sqr, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_ALR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
