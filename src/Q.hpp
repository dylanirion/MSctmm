#ifndef _Q_
#define _Q_

Rcpp::NumericMatrix getQ(const int nbStates, arma::vec rateparam, const time_t time, const double lng, const double lat, const int group, const Rcpp::String model);

#endif
