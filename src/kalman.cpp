#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "mat.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//' Kalman filter
//'
//' This code is adapted from the package ctmm (Calabrese et al., 2016) crawl (Johnson et al., 2008),
//' and MScrawl (Michelot and Blackwell, 2019).
//'
//' @name kalman_rcpp
//' @param data Matrix of data, including columns \code{"x"}, \code{"y"},
//' \code{"time"}, \code{"ID"} and \code{"state"} (in that order).
//' @param param Vector of movement parameters (\code{"tau_vel"}, \code{"tau_pos"}, and \code{"sigma"})
//' @param fixmu Vector of mean locations for the OUF process (\code{"x"}, \code{"y"})
//' @param Hmat Matrix of observation error variance (four columns, and one row
//' for each row of data)
//'
//' @return Log-likelihood
//'
//' @references
//' Calabrese, J.M., Fleming, C.H. and Gurarie, E. (2016).
//' ctmm: an r package for analyzing animal relocation data as a continuous‐time stochastic process.
//' Methods Ecol Evol, 7: 1124-1132. doi:10.1111/2041-210X.12559
//'
//' Fleming, C.H., Sheldon, D., Gurarie, E., Fagan, W.F., LaPoint, S., Calabrese, J.M. (2017).
//' Kálmán filters for continuous-time movement models.
//' Ecol Inform, 40: 8-21. doi:10.1016/j.ecoinf.2017.04.008
//'
//' Johnson, D.S., London, J.M., Lea, M.A., and Durban, J.W. (2008).
//' Continuous-time correlated random walk model for animal telemetry data.
//' Ecology, 89: 1208-1215. doi:10.1890/07-1032.1
//'
//' Michelot, T., Blackwell, P.G. (2019).
//' State‐switching continuous‐time correlated random walks.
//' Methods Ecol Evol, 10: 637-649. doi:10.1111/2041-210X.13154
//'
//' @export
// [[Rcpp::export]]
List kalman_rcpp( arma::mat& data, arma::vec param, arma::vec fixmu, arma::mat& Hmat ) {

  int nbData = data.n_rows;
  int N = nbData;
  int nbState = param.size() / 3;
  int nbID = 0;

  // unpack data
  mat X = data.cols( 0, 1 );    // x, y
  vec time = data.col( 2 );     // time
  vec ID = data.col( 3 );       // ID
  vec S = data.col( 4 );        // state
  vec dt( nbData, fill::ones ); // time intervals
  dt(0) = datum::inf;
  dt.subvec( 1, nbData - 1 ) = diff( time );

  mat aest( 4, 3, fill::zeros );  // state estimate
  //  x    mu_x    0
  //  vx   mu_vx   0
  //  y    0       mu_y
  //  vy   0       mu_vy
  mat Pest( 4, 4, fill::zeros );  //covariance estimate

  // unpack parameters
  vec tau_pos = param.subvec( 0, nbState - 1 );
  vec tau_vel = param.subvec( nbState, 2 * nbState - 1 );
  vec sigma = param.subvec( 2 * nbState, 3 * nbState - 1 );

  // define all empty matrices and vectors needed for the Kalman Filter and likelihood calculation
  mat Z{{1, 0 ,0 ,0}, {0, 0, 1, 0}};   // observation model which maps the true state space into the observed space (P, Hk)
  mat I(4, 4, fill::eye);              // identity matrix
  mat H(2, 2, fill::zeros);            // the covariance of the observation noise (error, Rk)
  mat T(4, 4);                         // state transition model which is applied to the previous state xk−1 (Green, Fk)
  mat Q(4, 4);                         // covariance of the process noise (Sigma, Q)
  mat K(4, 2, fill::zeros);            // Kalman Gain ( Gain, Kk )
  cube u(2, 3, nbData, fill::zeros);   // measurement residual (zRes, yk)
  cube iF(2, 2, nbData, fill::zeros);  // inverse of residual covariance, residual "precision" (isRes)
  cube uiF(2, 2, nbData, fill::zeros); // measurement residual * inverse of residual covariance (uisRes)
  cube zRes(2, 1, nbData, fill::zeros);// final measurement residual
  cube ziF(2, 1, nbData, fill::zeros); // final measurement residual * inverse of residual covariance
  cube mu(2, 1, nbData, fill::zeros);
  colvec logdetF(nbData, fill::zeros); // log determinant of residual covariance
  mat mu_out(nbState, 2);
  mu_out.fill(NA_REAL);

  // vector to keep track of the indices for first positions that are IOU, and size
  uvec iou(nbData);
  int k = 0;

  double logDet;
  double llk;

  // Kalman filter iterations
  for( int i = 0; i < nbData; i++ ) {

    // if first location or new individual
    if( i == 0 || ID( i ) != ID( i - 1 ) ) {
      nbID++;
      // if starting in IOU, keep track of index
      if( std::isinf( tau_pos( S(i) - 1 ) ) ) {
        iou(k) = i;
        k++;
      }
      // reset dt to inf (when new individual)
      dt(i) = datum::inf;

      // initialise state mean
      aest = zeros(4, 3);
      // and initial state covariance matrix
      Pest = makeQ(tau_pos(S(i) - 1), tau_vel(S(i) - 1), sigma(S(i) - 1), dt(i));
    }

    // update our estimate (if, missing obs skip to prediction)
    if(!R_IsNA(X(i, 0))) {
      H(0, 0) = Hmat(i, 0);
      H(1, 1) = Hmat(i, 1);
      H(0, 1) = Hmat(i, 2);
      H(1, 0) = Hmat(i, 3);
      mat aobs = join_rows(X.row(i).t(), eye(2, 2));
      // measurement residual (zRes, u)
      u.slice(i) = aobs - (Z * aest);

      // residual covariance (sRes, F)
      mat PestZt = Pest * Z.t();
      PestZt = PestZt.replace(datum::nan, 0);
      mat ZPestZt = Z * PestZt;
      ZPestZt = ZPestZt.replace(datum::nan, 0);
      mat F = ZPestZt + H;                                //residual covariance (sRes, Sk)
      iF.slice(i) = F.i();
      logdetF(i) = det(F) > 0 ? log(std::abs(det(F))) : datum::inf;
      uiF.slice(i) = iF.slice(i) * u.slice(i).submat(0, 1, 1, 2);
      uiF.slice(i) = uiF.slice(i).replace(datum::nan, 0);

      // Kalman gain
      K = PestZt * iF.slice(i);
      //if gain is inf or nan, replace with Z
      uvec idx = find_nonfinite(K);
      mat Zt = Z.t();
      K(idx) = Zt(idx);

      // update concurrent state estimate (zCon, aest)
      aest = aest + ( K * u.slice(i) );
      // update concurrent covariance estimate (sCon, Pest)
      mat J = I - ( K * Z );
      mat JPest = J * Pest;
      JPest = JPest.replace( datum::nan, 0 );
      mat KHKt = K * H * K.t();
      KHKt = KHKt.replace( datum::nan, 0 );
      Pest = JPest * J.t() + KHKt;
    }


    // proceed with forecast
    if( i < nbData - 1 && ID( i ) == ID( i + 1 ) ) {
      T = makeT( tau_pos( S( i + 1 ) - 1 ), tau_vel( S( i + 1 ) - 1), dt( i + 1 ) );
      Q = makeQ( tau_pos( S( i + 1 ) - 1 ), tau_vel( S( i + 1 ) - 1), sigma( S( i + 1 ) - 1 ), dt( i + 1 ) );

      //predict state estimate (zFor, aest)
      aest = T * aest;

      //predict covariance estimate (sFor, Pest)
      mat tcon = Pest;
      bool any = tcon.has_inf();
      if( any ) {
        tcon = tcon.replace( datum::inf, 0 );
      }

      Pest = T * tcon * T.t() + Q;

      if( any ) {
        vec Pdiag = Pest.diag();
        bool anyP = Pdiag.has_inf();
        uvec idx(4);
        int l = 0;
        if( anyP ) {
          for( int j = 0; j < Pdiag.size(); j++ ){
            if( std::isinf( Pdiag(j) ) ){
              idx(l) = j;
              l++;
            }
          }
          Pest.rows( idx.head(l) ).zeros();
          Pest.cols( idx.head(l) ).zeros();
          vec Pdiag2 = Pest.diag();
          Pdiag2( idx.head(l) ).fill( datum::inf );
          Pest.diag() = Pdiag2;
        }
      }
    }
  } // end filter iterations

  // shed all missing (skipped) obs
  uvec na_xy = find_nonfinite(X.col(0));
  u.shed_slices(na_xy);
  iF.shed_slices(na_xy);
  uiF.shed_slices(na_xy);
  ziF.shed_slices(na_xy);
  zRes.shed_slices(na_xy);
  mu.shed_slices(na_xy);
  X.shed_rows(na_xy);
  S.shed_rows(na_xy);

  // drop logDet for missing data and starting IOUs
  uvec drop = join_cols(na_xy, iou.head(k));
  logdetF.shed_rows(drop);
  N -= drop.size();

  // calculate state based mu
  // @todo: is element wise multiply on cube slower than reshaping? is there some kind of tensor operation (contraction?)
  for(int i = 0; i < nbState; i++) {

    // if IOU state, use first location for mu
    // @todo: do I need to do this per 'bout'
    if(std::isinf(tau_pos(i))){
      for(int j = 0; j < nbID; j++) {
        uvec iou_idx = find(S == i + 1);
        uvec ind_idx = find(ID == j + 1);
        uvec idx = intersect(iou_idx, ind_idx);
        if(idx.size() >= 1) {
          mu.each_slice(idx) = X.row(idx.min()).t();
        }
      }


    // otherwise calculate mu from residuals (maximum likelihood solution)
    } else {
      uvec ouf_idx = find(S == i + 1);        // Find indices where state i occurs
      if(R_IsNA(fixmu(i * 2)) && R_IsNA(fixmu(i * 2 + 1))) {
        cube uiF_ouf = uiF.slices(ouf_idx);   // measurement residual * inverse of residual covariance
        cube u_ouf = u.slices(ouf_idx);       // measurement residual
        cube D(2, 1, 1);                      // 'B' in eqns B30, B31
        D.tube(0,0) = sum(uiF_ouf.tube(0,0) % u_ouf.tube(0,0), 2); // data x
        D.tube(1,0) = sum(uiF_ouf.tube(1,1) % u_ouf.tube(1,0), 2); // data y
        cube W = sum(uiF_ouf % u_ouf.tube(0,1,1,2), 2);            // mean x, y
        //mat mu_m = inv(W.slice(0)) * D.slice(0);
        mat mu_m = solve(W.slice(0), eye(2,2)) * D.slice(0);
        mu.each_slice(ouf_idx) = mu_m;
        //save mu to output
        mu_out(i, 0) = mu_m(0,0);
        mu_out(i, 1) = mu_m(1,0);
      } else {
        mu.each_slice(ouf_idx) = fixmu.subvec(i * 2, i * 2 + 1);
        //save mu to output
        mu_out(i, 0) = fixmu(i * 2);
        mu_out(i, 1) = fixmu(i * 2 + 1);
      }
    }
  } // end mu

  // detrend
  for( int i = 0; i < zRes.n_slices; ++i ) {
    zRes.slice(i) = u.slice(i).col(0) - u.slice(i).cols(1,2) * mu.slice(i);
    ziF.slice(i) = iF.slice(i) * zRes.slice(i);
    ziF.slice(i) = ziF.slice(i).replace( datum::nan, 0 );
  } //end detrend

  double sigmaK = as_scalar( sum( sum( ziF % zRes, 2 ) ) ) / ( 2 * N );

  logDet = mean( logdetF );

  llk = -1 * ( sigmaK - 1 ) - logDet / 2;
  llk = N * ( llk + ( ( -1 * log( 2 * M_PI ) - 1 ) ) );
  llk = std::isnan(llk) ? -1 * datum::inf : llk;

  return List::create(
    Rcpp::Named("llk") = llk,
    Rcpp::Named("mu") = mu_out
  );
}
