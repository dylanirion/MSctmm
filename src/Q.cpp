#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Make Q matrix
//
Rcpp::NumericMatrix getQ(const int nbStates, const double time) {

  mat Q(nbStates, nbStates, fill::ones);
  Q.diag() = Q.diag() * -1;

  return Rcpp::wrap(Q);
}
