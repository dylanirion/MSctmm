#ifndef _SMOOTH_
#define _SMOOTH_
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "kalman.hpp"
#include "pdsolve.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List smooth_rcpp(const arma::mat &data, const int &nbStates, const arma::vec &param, const arma::vec &fixmu, const arma::mat &Hmat);

#endif
