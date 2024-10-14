#ifndef _MAT_
#define _MAT_
#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat makeT(double tau_pos, double tau_vel, double dt);
arma::mat makeQ(double tau_pos, double tau_vel, arma::mat sigma, double dt);
Rcpp::NumericVector dexp2(Rcpp::NumericVector x, Rcpp::Nullable<Rcpp::NumericVector> Exp = R_NilValue);
Rcpp::NumericVector dexp1(Rcpp::NumericVector x, Rcpp::Nullable<Rcpp::NumericVector> Exp = R_NilValue);
Rcpp::NumericVector sinc(Rcpp::NumericVector x, Rcpp::Nullable<Rcpp::NumericVector> SIN = R_NilValue);
Rcpp::NumericVector sinch(Rcpp::NumericVector x, Rcpp::Nullable<Rcpp::NumericVector> SINH = R_NilValue);

#endif
