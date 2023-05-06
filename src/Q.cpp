#include <RcppArmadillo.h>
#include <ctime>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Make Q matrix
// [[Rcpp::export]]
Rcpp::NumericMatrix getQ(const int nbStates, arma::vec alpha, arma::vec t_alpha, const time_t time, const double lng, const double lat, const String model) {
  //TODO convert time to day of year
  tm *t = localtime(&time);
  int yday = t->tm_yday;

  mat Q(nbStates, nbStates, fill::ones);
  Q.diag() = Q.diag() * -1;

  if (model == "time_out_time_in") {
    // time-varying rate in and out
    //FB -> trans
    Q(0,4) = alpha(0)/(1+exp(-alpha(0) * (yday - t_alpha(0))));
    Q(0,0) = Q(0,4) * -1;
    //GB -> trans
    Q(1,4) = alpha(0)/(1+exp(-alpha(0) * (yday - t_alpha(0))));
    Q(1,1) = Q(1,4) * -1;
    //MB -> trans
    Q(2,4) = alpha(0)/(1+exp(-alpha(0) * (yday - t_alpha(0))));
    Q(2,2) = Q(2,4) * -1;
    //AB -> trans
    Q(3,4) = alpha(0)/(1+exp(-alpha(0) * (yday - t_alpha(0))));
    Q(3,3) = Q(3,4) * -1;
    // trans -> FB
    Q(4,0) = alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1))))/4;
    // trans -> GB
    Q(4,1) = alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1))))/4;
    // trans -> MB
    Q(4,2) = alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1))))/4;
    // trans -> AB
    Q(4,3) = alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1))))/4;
    Q(4,4) = -(alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1)))));
    //Impossible transitions
    Q(0,1) = 0; // This might actually be possible (GB->FB)
    Q(0,2) = 0;
    Q(0,3) = 0;
    Q(1,0) = 0; // This might actually be possible (FB->GB)
    Q(1,2) = 0;
    Q(1,3) = 0;
    Q(2,0) = 0;
    Q(2,1) = 0;
    Q(2,3) = 0;
    Q(3,0) = 0;
    Q(3,1) = 0;
    Q(3,2) = 0;
  }

  return Rcpp::wrap(Q);
}
