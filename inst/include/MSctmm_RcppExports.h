// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_MSctmm_RCPPEXPORTS_H_GEN_
#define RCPP_MSctmm_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace MSctmm {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("MSctmm", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("MSctmm", "_MSctmm_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in MSctmm");
            }
        }
    }

    inline arma::mat sample_path_mr_(const int a, const int b, const double t0, const double t1, const double lng0, const double lat0, const double lng1, const double lat1, const int group, const double k, const int& nbStates, const arma::vec& param, const arma::vec& mu, const arma::mat& Hmat, Rcpp::List rateparam, const String model) {
        typedef SEXP(*Ptr_sample_path_mr_)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_path_mr_ p_sample_path_mr_ = NULL;
        if (p_sample_path_mr_ == NULL) {
            validateSignature("arma::mat(*sample_path_mr_)(const int,const int,const double,const double,const double,const double,const double,const double,const int,const double,const int&,const arma::vec&,const arma::vec&,const arma::mat&,Rcpp::List,const String)");
            p_sample_path_mr_ = (Ptr_sample_path_mr_)R_GetCCallable("MSctmm", "_MSctmm_sample_path_mr_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_path_mr_(Shield<SEXP>(Rcpp::wrap(a)), Shield<SEXP>(Rcpp::wrap(b)), Shield<SEXP>(Rcpp::wrap(t0)), Shield<SEXP>(Rcpp::wrap(t1)), Shield<SEXP>(Rcpp::wrap(lng0)), Shield<SEXP>(Rcpp::wrap(lat0)), Shield<SEXP>(Rcpp::wrap(lng1)), Shield<SEXP>(Rcpp::wrap(lat1)), Shield<SEXP>(Rcpp::wrap(group)), Shield<SEXP>(Rcpp::wrap(k)), Shield<SEXP>(Rcpp::wrap(nbStates)), Shield<SEXP>(Rcpp::wrap(param)), Shield<SEXP>(Rcpp::wrap(mu)), Shield<SEXP>(Rcpp::wrap(Hmat)), Shield<SEXP>(Rcpp::wrap(rateparam)), Shield<SEXP>(Rcpp::wrap(model)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

}

#endif // RCPP_MSctmm_RCPPEXPORTS_H_GEN_
