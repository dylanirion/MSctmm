// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/MSctmm.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getQ
Rcpp::NumericMatrix getQ(const int nbStates, arma::vec rateparam, const time_t time, const double lng, const double lat, const int group, const String model);
RcppExport SEXP _MSctmm_getQ(SEXP nbStatesSEXP, SEXP rateparamSEXP, SEXP timeSEXP, SEXP lngSEXP, SEXP latSEXP, SEXP groupSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nbStates(nbStatesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rateparam(rateparamSEXP);
    Rcpp::traits::input_parameter< const time_t >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double >::type lng(lngSEXP);
    Rcpp::traits::input_parameter< const double >::type lat(latSEXP);
    Rcpp::traits::input_parameter< const int >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const String >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(getQ(nbStates, rateparam, time, lng, lat, group, model));
    return rcpp_result_gen;
END_RCPP
}
// kalman_rcpp
List kalman_rcpp(const arma::mat& data, const int& nbStates, const arma::vec& param, const arma::vec& fixmu, const arma::mat& Hmat);
RcppExport SEXP _MSctmm_kalman_rcpp(SEXP dataSEXP, SEXP nbStatesSEXP, SEXP paramSEXP, SEXP fixmuSEXP, SEXP HmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type nbStates(nbStatesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type fixmu(fixmuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_rcpp(data, nbStates, param, fixmu, Hmat));
    return rcpp_result_gen;
END_RCPP
}
// sample_path_mr_
arma::mat sample_path_mr_(const int a, const int b, const double t0, const double t1, const double lng0, const double lat0, const double lng1, const double lat1, const int group, const double k, const int& nbStates, const arma::vec& param, const arma::vec& mu, const arma::mat& Hmat, Rcpp::List rateparam, const String model);
static SEXP _MSctmm_sample_path_mr__try(SEXP aSEXP, SEXP bSEXP, SEXP t0SEXP, SEXP t1SEXP, SEXP lng0SEXP, SEXP lat0SEXP, SEXP lng1SEXP, SEXP lat1SEXP, SEXP groupSEXP, SEXP kSEXP, SEXP nbStatesSEXP, SEXP paramSEXP, SEXP muSEXP, SEXP HmatSEXP, SEXP rateparamSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int >::type a(aSEXP);
    Rcpp::traits::input_parameter< const int >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< const double >::type lng0(lng0SEXP);
    Rcpp::traits::input_parameter< const double >::type lat0(lat0SEXP);
    Rcpp::traits::input_parameter< const double >::type lng1(lng1SEXP);
    Rcpp::traits::input_parameter< const double >::type lat1(lat1SEXP);
    Rcpp::traits::input_parameter< const int >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const double >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int& >::type nbStates(nbStatesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type rateparam(rateparamSEXP);
    Rcpp::traits::input_parameter< const String >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_path_mr_(a, b, t0, t1, lng0, lat0, lng1, lat1, group, k, nbStates, param, mu, Hmat, rateparam, model));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _MSctmm_sample_path_mr_(SEXP aSEXP, SEXP bSEXP, SEXP t0SEXP, SEXP t1SEXP, SEXP lng0SEXP, SEXP lat0SEXP, SEXP lng1SEXP, SEXP lat1SEXP, SEXP groupSEXP, SEXP kSEXP, SEXP nbStatesSEXP, SEXP paramSEXP, SEXP muSEXP, SEXP HmatSEXP, SEXP rateparamSEXP, SEXP modelSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_MSctmm_sample_path_mr__try(aSEXP, bSEXP, t0SEXP, t1SEXP, lng0SEXP, lat0SEXP, lng1SEXP, lat1SEXP, groupSEXP, kSEXP, nbStatesSEXP, paramSEXP, muSEXP, HmatSEXP, rateparamSEXP, modelSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// smooth_rcpp
List smooth_rcpp(const arma::mat& data, const int& nbStates, const arma::vec& param, const arma::vec& fixmu, const arma::mat& Hmat);
RcppExport SEXP _MSctmm_smooth_rcpp(SEXP dataSEXP, SEXP nbStatesSEXP, SEXP paramSEXP, SEXP fixmuSEXP, SEXP HmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type nbStates(nbStatesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type fixmu(fixmuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_rcpp(data, nbStates, param, fixmu, Hmat));
    return rcpp_result_gen;
END_RCPP
}
// makeMu
arma::mat makeMu(double tau_pos, double tau_vel, double dt);
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
arma::mat makeSigma(double tau_pos, double tau_vel, arma::mat sigma, double dt);
RcppExport SEXP _MSctmm_makeSigma(SEXP tau_posSEXP, SEXP tau_velSEXP, SEXP sigmaSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tau_pos(tau_posSEXP);
    Rcpp::traits::input_parameter< double >::type tau_vel(tau_velSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(makeSigma(tau_pos, tau_vel, sigma, dt));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _MSctmm_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::mat(*sample_path_mr_)(const int,const int,const double,const double,const double,const double,const double,const double,const int,const double,const int&,const arma::vec&,const arma::vec&,const arma::mat&,Rcpp::List,const String)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _MSctmm_RcppExport_registerCCallable() { 
    R_RegisterCCallable("MSctmm", "_MSctmm_sample_path_mr_", (DL_FUNC)_MSctmm_sample_path_mr__try);
    R_RegisterCCallable("MSctmm", "_MSctmm_RcppExport_validate", (DL_FUNC)_MSctmm_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_MSctmm_getQ", (DL_FUNC) &_MSctmm_getQ, 7},
    {"_MSctmm_kalman_rcpp", (DL_FUNC) &_MSctmm_kalman_rcpp, 5},
    {"_MSctmm_sample_path_mr_", (DL_FUNC) &_MSctmm_sample_path_mr_, 16},
    {"_MSctmm_smooth_rcpp", (DL_FUNC) &_MSctmm_smooth_rcpp, 5},
    {"_MSctmm_makeMu", (DL_FUNC) &_MSctmm_makeMu, 3},
    {"_MSctmm_makeSigma", (DL_FUNC) &_MSctmm_makeSigma, 4},
    {"_MSctmm_RcppExport_registerCCallable", (DL_FUNC) &_MSctmm_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_MSctmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
