
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "mat.hpp"

//' Kalman filter
//'
//' This code is adapted from the package ctmm (Calabrese et al., 2016) crawl (Johnson et al., 2008),
//' and MScrawl (Michelot and Blackwell, 2019).
//'
//' @name kalman_rcpp
//' @param data Matrix of data, including columns \code{"x"}, \code{"y"},
//' \code{"time"}, \code{"ID"} and \code{"state"} (in that order).
//' @param param Vector of movement parameters (\code{"tau_vel"}, \code{"tau_pos"}, and \code{"sigma"})
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
Rcpp::NumericVector kalman_rcpp( arma::mat& data, arma::vec param, arma::mat& Hmat ) {

  int nbData = data.n_rows;
  int N = nbData;
  int nbState = param.size() / 3;

  // unpack data
  mat X = data.cols( 0, 1 );  // x, y
  vec time = data.col( 2 );   // time
  vec ID = data.col( 3 );     // ID
  vec S = data.col( 4 );      // state
  // time intervals
  vec dt( nbData, fill::ones );
  dt(0) = datum::inf;
  dt.subvec( 1, nbData - 1 ) = diff( time );
  mat aest( 4, 3, fill::zeros );  // state estimate
  mat Pest( 4, 4, fill::zeros );  //covariance estimate

  // unpack parameters
  vec tau_pos = param.subvec( 0, nbState - 1 );
  vec tau_vel = param.subvec( nbState, 2 * nbState - 1 );
  vec sigma = param.subvec( 2 * nbState, 3 * nbState - 1 );

  // define all matrices and vectors needed below
  mat Z( 2, 4, fill::zeros ); // observation model which maps the true state space into the observed space ( P, Hk )
  Z( 0, 0 ) = 1;
  Z( 1, 2 ) = 1;
  mat I( 4, 4, fill::eye );   // identity matrix
  mat H( 2, 2, fill::zeros ); // the covariance of the observation noise ( error, Rk )
  mat T( 4, 4 );              // state transition model which is applied to the previous state xk−1 ( Green, Fk )
  mat Q( 4, 4 );              // covariance of the process noise ( Sigma, Q )
  mat K( 4, 2, fill::zeros ); // Kalman Gain ( Gain, Kk )
  cube u( 2, 3, nbData, fill::zeros );  // measurement residual ( zRes, yk )
  mat zRes( 2 * nbData, 1, fill::zeros ); // final measurement residual
  cube iF( 2, 2, nbData, fill::zeros ); // inverse of residual covariance ( isRes )
  cube uiF( 2, 2, nbData, fill::zeros ); // measurement residual * inverse of residual covariance ( uisRes )
  colvec logdetF( nbData, fill::zeros ); // log determinant of residual covariance
  mat mu( 2, nbState);
  mu.fill(datum::inf);
  mat iou( 1, 2, fill::zeros );  //matrix to track first and last index of each IOU bout

  double logDet;
  double llk;

  // Kalman filter iterations
  for( int i = 0; i < nbData; i++ ) {

    if( i == 0 || ID( i ) != ID( i - 1 ) ) {
      dt(i) = datum::inf;
      // If first location of track, initialise state mean
      aest = zeros( 4, 3 );
      // and initial state covariance matrix
      Pest = zeros( 4, 4 );
      Pest = makeQ( tau_pos( S(i) - 1 ), tau_vel( S(i) - 1 ), sigma( S(i) - 1 ), dt(i) );
    }

    // if IOU, and first position of new bout or new individual, track start index and reduce N by 1
    if( std::isinf( tau_pos( S(i) - 1 ) ) ) {
      if( i == 0 || S(i) != S( i - 1 ) || ID( i ) != ID( i - 1 ) ) {
        N = N - 1;
        iou( iou.n_rows - 1, 0 ) = i;
      } else if( i == nbData - 1 ) { //if still IOU, but last iteration, track end index
        iou( iou.n_rows - 1, 1 ) = i - 1;
      }
    } else if( !std::isinf( tau_pos( S(i) - 1 ) ) && i > 0 ) {
      if( S(i) != S( i - 1 ) || ID( i ) != ID( i - 1 ) )  {
        iou( iou.n_rows - 1, 1 ) = i - 1;
        if( i != nbData - 1) {
          iou.resize( iou.n_rows + 1, 2 );
        }
      }
    }

    if( i < nbData - 1 ) {
      T = makeT( tau_pos( S( i + 1 ) - 1 ), tau_vel( S( i + 1 ) - 1), dt( i + 1 ) );
      Q = makeQ( tau_pos( S( i + 1 ) - 1 ), tau_vel( S( i + 1 ) - 1), sigma( S( i + 1 ) - 1 ), dt( i + 1 ) );
    }

    if( R_IsNA( X( i,0 ) ) ) {
      // if missing observation, we just update/forecast
      aest = T * aest;
      Pest = T * Pest * T.t() + Q;
    } else {
      H( 0, 0 ) = Hmat( i, 0 );
      H( 1, 1 ) = Hmat( i, 1 );
      H( 0, 1 ) = Hmat( i, 2 );
      H( 1, 0 ) = Hmat( i, 3 );
      mat uslice( 2, 3, fill::zeros );
      uslice.col(0) = X.row(i).t();
      uslice( 0, 1 ) = 1;
      uslice( 1, 2 ) = 1;
      // measurement residual (zRes, u)
      u.slice(i) = uslice - ( Z * aest );                  //(2,3)
      // residual covariance (sRes, F)
      mat PestZt = Pest * Z.t();                          //(4,2)
      PestZt = PestZt.replace( datum::nan, 0 );
      mat ZPestZt = Z * PestZt;                           //(2,2)
      ZPestZt = ZPestZt.replace( datum::nan, 0 );
      mat F = ZPestZt + H;                                // (2,2) residual covariance ( sRes, Sk )
      iF.slice(i) = F.i();
      logdetF(i) = det(F) > 0 ? log( std::abs( det(F) ) ) : datum::inf;
      uiF.slice(i) =  iF.slice(i) * u.slice(i).submat( 0, 1, 1, 2 );    //(2,2)
      uiF.slice(i) = uiF.slice(i).replace( datum::nan, 0 );

      // Kalman gain
      K = PestZt * iF.slice(i);                                    // (4,2)
      //if gain is inf or nan, replace with Z
      uvec idx = find_nonfinite( K );
      mat Zt = Z.t();
      K(idx) = Zt(idx);

      // concurrent state estimate (zCon, aest)
      aest = aest + ( K * u.slice(i) );
      // concurrent covariance estimate (sCon, Pest)
      mat J = I - ( K * Z );              //(4,4)
      mat JPest = J * Pest;               //(4,4)
      JPest = JPest.replace( datum::nan, 0 );
      mat KHKt = K * H * K.t();           //(4,4)
      KHKt = KHKt.replace( datum::nan, 0 );
      Pest = JPest * J.t() + KHKt;

      //update state estimate (zFor, aest)
      aest = T * aest;

      //update covariance estimate (sFor, Pest)
      mat tcon = Pest;
      bool any = tcon.has_inf();
      if(any) {
        tcon = tcon.replace( datum::inf, 0 );
      }

      Pest = T * tcon * T.t() + Q;

      if(any) {
        vec Pdiag = Pest.diag();
        bool anyP = Pdiag.has_inf();
        uvec idx;
        if(anyP) {
          for(int j=0; j<Pdiag.size(); j++){
            if( std::isinf( Pdiag(j) ) ){
              int sz = idx.size();
              idx.resize( sz + 1 );
              idx(sz) = j;
            }
          }
          Pest.rows( idx ).zeros();
          Pest.cols( idx ).zeros();
          vec Pdiag2 = Pest.diag();
          Pdiag2( idx ).fill( datum::inf );
          Pest.diag() = Pdiag2;
        }
      }
    }
  }

  // u and uiF have dimensions (2, 3, nbData) and (2, 2, nbData)
  // then get permuted to (nbData, 2, 3) and (2, nbData, 2)
  // and reshaped to (2x nbData, 3, 1) and (2x nbData, 2, 1)

  //permute u(r,c,s) to u(s,r,c) and reshape
  cube uperm( nbData, 2, 3 );
  for( int c = 0; c < u.n_cols; ++c )
    for( int r = 0; r < u.n_rows; ++r )
      for( int s = 0; s < u.n_slices; ++s )
        uperm( s, r, c ) = u( r, c, s );
  uperm.reshape( 2 * nbData, 3, 1 );

  //permute uiF(r,c,s) to uiF(c,s,r) and reshape
  cube uiFperm( 2, nbData, 2 );
  for (int c = 0; c < uiF.n_cols; ++c)
    for (int r = 0; r < uiF.n_rows; ++r)
      for (int s = 0; s < uiF.n_slices; ++s)
        uiFperm(c, s, r) = uiF(r, c, s);
  uiFperm.reshape(2, nbData * 2, 1 );

  // calculate state based mu
  for( int i = 0; i < nbState; i++ ) {
    mat uiFslice = uiFperm.slice(0);
    mat uSlice = uperm.slice(0);
    uvec col( 1, fill::zeros );
    uvec cols = regspace<uvec>(1,2);

    // if IOU overwrite NaN stationary mean value with first location of each IOU bout
    if( std::isinf( tau_pos( i ) ) ){
      for(int j = 0; j < iou.n_rows; j++ ){
        colvec iou_mu = u.slice( iou(j,0) ).col(0);
        uvec iou_ids = sort( join_cols ( regspace<uvec>( iou(j,0), iou(j,1) ), regspace<uvec>( iou(j,0) + nbData, iou(j,1) + nbData ) ) );

        zRes.rows(iou_ids) = uSlice.submat( iou_ids, col ) - ( uSlice.submat( iou_ids, cols ) * iou_mu );
      }
    } else {
      uvec ids = find( S == i + 1 );        // Find indices where state i occurs
      ids = sort( join_cols ( ids, ids + nbData ) );
      mat D = uiFslice.cols(ids) * uSlice.submat( ids, col );    //(2,1)
      mat W = uiFslice.cols(ids) * uSlice.submat( ids, cols );   //(2,2)
      mat iW( 2, 2, fill::zeros );
      iW = W.i();
      mu.col(i) = iW * D;  // (2,1)

      zRes.rows(ids) = uSlice.submat( ids, col ) - ( uSlice.submat( ids, cols ) * mu.col(i) );
    }
  }

  zRes.reshape(nbData,2);
  //mu.reshape(nbState,2);
  cube ziF(2,1,nbData, fill::zeros);
  for (int s = 0; s < iF.n_slices; ++s) {
    ziF.slice(s) = iF.slice(s) * zRes.row(s).t();       //(2,1) Line 450
  }
  //permute ziF(r,c,s) to ziF(c,s,r)
  cube ziFperm(1,nbData,2);
  for (int c = 0; c < ziF.n_cols; ++c)
    for (int r = 0; r < ziF.n_rows; ++r)
      for (int s = 0; s < ziF.n_slices; ++s)
        ziFperm(c, s, r) = ziF(r, c, s);
  ziFperm.reshape(1, nbData * 2, 1 );
  zRes.reshape( 2 * nbData, 1 );

  mat sigmaK = ( ziFperm.slice(0).replace( datum::nan, 0 )  * zRes ) / ( 2 * N );

  // if IOU, drop first logDet for each IOU bout because couldn't condition off initial state and will be Inf
  uvec drop = as<uvec>( wrap( iou.col(0) ) );
  if( tau_pos.has_inf() ) {
    logdetF.shed_rows( drop );
  }

  logDet = mean( logdetF );

  llk = -1 * ( as_scalar( sigmaK ) - 1 ) - logDet / 2;
  llk = N * ( llk + ( ( -1 * log( 2 * M_PI ) - 1 ) / N ) );
  llk = std::isnan(llk) ? -1 * datum::inf : llk;

  NumericVector out = as<NumericVector>( wrap( join_cols( as<colvec>( wrap( llk ) ), vectorise( mu ) ) ) );
  out.attr("dim") = R_NilValue;

  return out;

}
