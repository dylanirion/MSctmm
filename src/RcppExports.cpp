// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// kalman_rcpp
NumericVector kalman_rcpp(arma::mat& data, arma::vec param, arma::vec fixmu, arma::mat& Hmat);
RcppExport SEXP _MSctmm_kalman_rcpp(SEXP dataSEXP, SEXP paramSEXP, SEXP fixmuSEXP, SEXP HmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fixmu(fixmuSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hmat(HmatSEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_rcpp(data, param, fixmu, Hmat));
    return rcpp_result_gen;
END_RCPP
}
// makeMu
NumericMatrix makeMu(double tau_pos, double tau_vel, double dt);
RcppExport SEXP _MSctmm_makeMu(SEXP tau_posSEXP, SEXP tau_velSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tau_pos(tau_posSEXP);
    Rcpp::traits::input_parameter< double >::type tau_vel(tau_velSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(makeMu(tau_pos, tau_vel, dt));
    return rcpp_result_gen;
END_RCPP
}
// makeSigma
NumericMatrix makeSigma(double tau_pos, double tau_vel, double sigma, double dt);
RcppExport SEXP _MSctmm_makeSigma(SEXP tau_posSEXP, SEXP tau_velSEXP, SEXP sigmaSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tau_pos(tau_posSEXP);
    Rcpp::traits::input_parameter< double >::type tau_vel(tau_velSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(makeSigma(tau_pos, tau_vel, sigma, dt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MSctmm_kalman_rcpp", (DL_FUNC) &_MSctmm_kalman_rcpp, 4},
    {"_MSctmm_makeMu", (DL_FUNC) &_MSctmm_makeMu, 3},
    {"_MSctmm_makeSigma", (DL_FUNC) &_MSctmm_makeSigma, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_MSctmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
