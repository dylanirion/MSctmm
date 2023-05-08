#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "mat.hpp"

//' Make mu matrix
//'
//' This code is adapted from the langevin function in the R package ctmm (Calabrese et al., 2016).
//' It is an R wrapper for the internal makeT function
//'
//' @name makeMu
//' @param tau_pos Parameter \eqn{tau[pos]} of the movement process
//' @param tau_vel Parameter \eqn{tau[vel]} of the movement process
//' @param dt Time interval
//'
//' @return Green's Function matrix
//'
//' @references
//' Calabrese, J.M., Fleming, C.H. and Gurarie, E. (2016).
//' ctmm: an r package for analyzing animal relocation data as a continuous‐time stochastic process.
//' Methods Ecol Evol, 7: 1124-1132. doi:10.1111/2041-210X.12559
//'
//' @export
// [[Rcpp::export]]
arma::mat makeMu( double tau_pos, double tau_vel, double dt ) {
  return makeT( tau_pos, tau_vel, dt );
}




//' Make covariance matrix
//'
//' This code is adapted from the langevin function in the R package ctmm (Calabrese et al., 2016).
//' It is an R wrapper for the internal makeQ function
//'
//' @name makeSigma
//' @param tau_pos Parameter \eqn{tau[pos]} of the movement process
//' @param tau_vel Parameter \eqn{tau[vel]} of the movement process
//' @param sigma Parameter \eqn{sigma} of the movement process
//' @param dt Time interval
//'
//' @return Covariance Matrix
//'
//' @references
//' Calabrese, J.M., Fleming, C.H. and Gurarie, E. (2016).
//' ctmm: an r package for analyzing animal relocation data as a continuous‐time stochastic process.
//' Methods Ecol Evol, 7: 1124-1132. doi:10.1111/2041-210X.12559
//'
//' This code is adapted from the langevin function in the R package ctmm (Calabrese et al., 2016).
//' @export
// [[Rcpp::export]]
arma::mat makeSigma( double tau_pos, double tau_vel, arma::mat sigma, double dt ) {
  return makeQ( tau_pos, tau_vel, sigma, dt );
}
