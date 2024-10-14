#ifndef _KALMAN_
#define _KALMAN_
#include <RcppArmadillo.h>
#include <math.h>
#include "pdsolve.hpp"
#include "mat.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

template <typename AestType, typename PestType>
void initialize_state(
    const int &i, const arma::vec &S, const arma::vec &tau_pos, const arma::vec &tau_vel,
    const arma::cube &sigma, arma::vec &dt, AestType &&aest, PestType &&Pest)
{
  // reset dt to inf
  dt(i) = arma::datum::inf;
  // initialise state mean
  aest.zeros();
  // and initial state covariance matrix
  Pest = makeQ(tau_pos(S(i) - 1), tau_vel(S(i) - 1), sigma.slice(S(i) - 1), dt(i));
}

template void initialize_state(
    const int &, const arma::vec &, const arma::vec &, const arma::vec &,
    const arma::cube &, arma::vec &, arma::mat &, arma::mat &);

template void initialize_state(
    const int &, const arma::vec &, const arma::vec &, const arma::vec &,
    const arma::cube &, arma::vec &, arma::subview_row<double> &, arma::Mat<double> &);

void prepare_sigma(const arma::vec &param, const int &nbStates, arma::cube &sigma);

bool is_observation_missing(const arma::mat &X, const int &i);

#endif