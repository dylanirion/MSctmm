#ifndef _SMOOTH_
#define _SMOOTH_

Rcpp::List smooth_rcpp(const arma::mat& data, int nbStates, const arma::vec param, const arma::vec fixmu, const arma::mat& Hmat);

#endif
