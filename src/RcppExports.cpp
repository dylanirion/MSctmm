// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// kalman_rcpp
double kalman_rcpp(arma::mat& data, arma::vec param, arma::mat& Hmat);
RcppExport SEXP _MSctmm_kalman_rcpp(SEXP dataSEXP, SEXP paramSEXP, SEXP HmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hmat(HmatSEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_rcpp(data, param, Hmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MSctmm_kalman_rcpp", (DL_FUNC) &_MSctmm_kalman_rcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MSctmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
