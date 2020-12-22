#define ARMA_DONT_PRINT_ERRORS
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
NumericVector kalman_rcpp( arma::mat& data, arma::vec param, arma::vec fixmu, arma::mat& Hmat ) {

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
  //  x    mu_x    0
  //  vx   mu_vx   0
  //  y    0       mu_y
  //  vy   0       mu_vy
  mat Pest( 4, 4, fill::zeros );  //covariance estimate

  // unpack parameters
  vec tau_pos = param.subvec( 0, nbState - 1 );
  vec tau_vel = param.subvec( nbState, 2 * nbState - 1 );
  vec sigma = param.subvec( 2 * nbState, 3 * nbState - 1 );

  // define all empty matrices and vectors needed for the Kalman Filter and elsewhere
  mat Z { { 1, 0 ,0 ,0 }, { 0, 0, 1, 0 } }; // observation model which maps the true state space into the observed space ( P, Hk )
  mat I( 4, 4, fill::eye );   // identity matrix
  mat H( 2, 2, fill::zeros ); // the covariance of the observation noise ( error, Rk )
  mat T( 4, 4 );              // state transition model which is applied to the previous state xk−1 ( Green, Fk )
  mat Q( 4, 4 );              // covariance of the process noise ( Sigma, Q )
  mat K( 4, 2, fill::zeros ); // Kalman Gain ( Gain, Kk )
  cube u( 2, 3, nbData, fill::zeros );      // measurement residual ( zRes, yk )
  cube iF( 2, 2, nbData, fill::zeros ); // inverse of residual covariance, residual "precision" ( isRes )
  cube uiF( 2, 2, nbData, fill::zeros ); // measurement residual * inverse of residual covariance ( uisRes )
  cube zRes( 2, 1, nbData, fill::zeros ); // final measurement residual
  cube ziF( 2, 1, nbData, fill::zeros); // final measurement residual * inverse of residual covariance
  cube mu( 2, 1, nbData, fill::zeros);
  mat iou( 1, 2, fill::zeros );  //matrix to track first and last index of each IOU bout
  bool iou_open = false;
  colvec logdetF( nbData, fill::zeros ); // log determinant of residual covariance
  colvec out( 2 * nbState + 1 );
  out.fill( NA_REAL );

  double logDet;
  double llk;

  // Kalman filter iterations
  for( int i = 0; i < nbData; i++ ) {

    if( i == 0 || ID( i ) != ID( i - 1 ) ) {
      dt(i) = datum::inf;
      // If first location of track, initialise state mean
      aest = zeros( 4, 3 );
      // and initial state covariance matrix
      Pest = makeQ( tau_pos( S(i) - 1 ), tau_vel( S(i) - 1 ), sigma( S(i) - 1 ), dt(i) );
    }

    // if current state is IOU
    if( std::isinf( tau_pos( S(i) - 1 ) ) ) {
      // and first position or first position of new bout or individual
      if( i == 0 || S(i) != S( i - 1 ) || ID( i ) != ID( i - 1 ) ) {
        //reduce N by 1
        N = N - 1;
        //check if we need to close a bout and add a new one
        if( iou_open ) {
          // track end index
          iou( iou.n_rows - 1, 1 ) = i - 1;
          // increase size
          iou.resize( iou.n_rows + 1, 2 );
          iou_open = false;
        }
        // if it's an augmented position ( i.e. NA ), use the next index as start
        if( R_IsNA( X( i, 0 ) ) ) {
          iou( iou.n_rows - 1, 0 ) = i + 1;
          iou_open = true;
        } else {
          // track start index
          iou( iou.n_rows - 1, 0 ) = i;
          // track iou_open
          iou_open = true;
        }
      // or last iteration, track end index
      } else if( i == nbData - 1 ) {
        iou( iou.n_rows - 1, 1 ) = i - 1;
        iou_open = false;
      }
    // else if not IOU
    } else if( !std::isinf( tau_pos( S(i) - 1 ) ) && i > 0 ) {
      // but the state switched or the individual changed
      if( S(i) != S( i - 1 ) || ID( i ) != ID( i - 1 ) )  {
        //close the last bout index
        iou( iou.n_rows - 1, 1 ) = i - 1;
        iou_open = false;
      }
    }

    if( i < nbData - 1 ) {
      T = makeT( tau_pos( S( i + 1 ) - 1 ), tau_vel( S( i + 1 ) - 1), dt( i + 1 ) );
      Q = makeQ( tau_pos( S( i + 1 ) - 1 ), tau_vel( S( i + 1 ) - 1), sigma( S( i + 1 ) - 1 ), dt( i + 1 ) );
    }

    // if missing observation, we can skip this and just update/forecast
    if( !R_IsNA( X( i, 0 ) ) ) {
      H( 0, 0 ) = Hmat( i, 0 );
      H( 1, 1 ) = Hmat( i, 1 );
      H( 0, 1 ) = Hmat( i, 2 );
      H( 1, 0 ) = Hmat( i, 3 );
      mat aobs = join_rows( X.row(i).t(), eye( 2, 2 ) );
      // measurement residual (zRes, u)
      u.slice(i) = aobs - ( Z * aest );

      // residual covariance (sRes, F)
      mat PestZt = Pest * Z.t();
      PestZt = PestZt.replace( datum::nan, 0 );
      mat ZPestZt = Z * PestZt;
      ZPestZt = ZPestZt.replace( datum::nan, 0 );
      mat F = ZPestZt + H;                                //residual covariance ( sRes, Sk )
      iF.slice(i) = F.i();
      logdetF(i) = det(F) > 0 ? log( std::abs( det(F) ) ) : datum::inf;
      uiF.slice(i) = iF.slice(i) * u.slice(i).submat( 0, 1, 1, 2 );
      uiF.slice(i) = uiF.slice(i).replace( datum::nan, 0 );

      // Kalman gain
      K = PestZt * iF.slice(i);
      //if gain is inf or nan, replace with Z
      uvec idx = find_nonfinite( K );
      mat Zt = Z.t();
      K(idx) = Zt(idx);

      // concurrent state estimate (zCon, aest)
      aest = aest + ( K * u.slice(i) );
      // concurrent covariance estimate (sCon, Pest)
      mat J = I - ( K * Z );
      mat JPest = J * Pest;
      JPest = JPest.replace( datum::nan, 0 );
      mat KHKt = K * H * K.t();
      KHKt = KHKt.replace( datum::nan, 0 );
      Pest = JPest * J.t() + KHKt;
    }

    //update state estimate (zFor, aest)
    aest = T * aest;

    //update covariance estimate (sFor, Pest)
    mat tcon = Pest;
    bool any = tcon.has_inf();
    if( any ) {
      tcon = tcon.replace( datum::inf, 0 );
    }

    Pest = T * tcon * T.t() + Q;

    if( any ) {
      vec Pdiag = Pest.diag();
      bool anyP = Pdiag.has_inf();
      uvec idx;
      if( anyP ) {
        for( int j = 0; j < Pdiag.size(); j++ ){
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

  // calculate state based mu
  // @todo: is element wise multiply on cube slower than reshaping? is there some kind of tensor operation (contraction?)
  for( int i = 0; i < nbState; i++ ) {

    // if IOU use first location of each IOU bout for mu
    if( std::isinf( tau_pos( i ) ) ){
      for(int j = 0; j < iou.n_rows; j++ ){
        uvec iou_ids = regspace<uvec>( iou(j,0), iou(j,1) );
        mu.each_slice(iou_ids) = X.row( iou(j,0) ).t();
      }
    // otherwise calculate mu from residuals
    } else {
      uvec ouf_ids = find( S == i + 1 );        // Find indices where state i occurs
      if( R_IsNA( fixmu( i * 2 ) ) && R_IsNA( fixmu( i * 2 + 1 ) ) ) {
        cube uiF_ouf = uiF.slices(ouf_ids);
        cube u_ouf = u.slices(ouf_ids);
        cube D( 2, 1, 1 );
        D.tube(0,0) = sum( uiF_ouf.tube(0,0) % u_ouf.tube(0,0), 2 );
        D.tube(1,0) = sum( uiF_ouf.tube(1,1) % u_ouf.tube(1,0), 2 );
        cube W = sum( uiF_ouf % u_ouf.tube(0,1,1,2), 2 );
        //mat mu_m = inv( W.slice(0) ) * D.slice(0);
        mat mu_m = solve( W.slice(0), eye(2,2) ) * D.slice(0);
        mu.each_slice(ouf_ids) = mu_m;
        out( i * 2 + 1 ) = mu_m(0,0);
        out( i * 2 + 2 ) =  mu_m(1,0);
      } else {
        mu.each_slice(ouf_ids) = fixmu.subvec( i * 2, i * 2 + 1 );
        out( i * 2 + 1 ) = fixmu( i * 2 );
        out( i * 2 + 2 ) = fixmu( i * 2 + 1 );
      }

    }
  }

  for( int i = 0; i < zRes.n_slices; ++i ) {
    zRes.slice(i) = u.slice(i).col(0) - u.slice(i).cols(1,2) * mu.slice(i);
    ziF.slice(i) = iF.slice(i) * zRes.slice(i);
    ziF.slice(i) = ziF.slice(i).replace( datum::nan, 0 );
  }
  uvec na_xy = find_nonfinite( X.col(0) );
  ziF.shed_slices( na_xy );
  zRes.shed_slices( na_xy );
  N -= na_xy.size();

  double sigmaK = as_scalar( sum( sum( ziF % zRes, 2 ) ) ) / ( 2 * N );

  // if IOU, drop first logDet for each IOU bout because couldn't condition off initial state and will be Inf
  if( tau_pos.has_inf() ) {
    uvec drop_iou = as<uvec>( wrap( iou.col(0) ) ); // is there a better way to do this?
    logdetF.shed_rows( drop_iou );
  }

  logDet = mean( logdetF );

  llk = -1 * ( sigmaK - 1 ) - logDet / 2;
  llk = N * ( llk + ( ( -1 * log( 2 * M_PI ) - 1 ) / N ) );
  llk = std::isnan(llk) ? -1 * datum::inf : llk;

  out(0) = llk;
  NumericVector output = as<NumericVector>( wrap( out ) );
  output.attr("dim") = R_NilValue;

  return output;

}
