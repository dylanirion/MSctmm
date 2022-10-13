#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Make Q matrix
//
Rcpp::NumericMatrix getQ(const int nbStates, const double time, const double long, const double lat, const std::string model) {

  mat Q(nbStates, nbStates, fill::ones);
  Q.diag() = Q.diag() * -1;

  if (model == "time_out_loc_in") {
    // time-varying rate out
    // location-varying rate in
    //FB -> trans
    Q(0,4) = 1;
    //GB -> trans
    Q(1,4) = 1;
    //MB -> trans
    Q(2,4) = 1;
    //AB -> trans
    Q(3,4) = 1;
  }

  return Rcpp::wrap(Q);
}
