#include <RcppArmadillo.h>
#include "mat.hpp"
#include "pdsolve.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' Kalman filter
//'
//' This code is adapted from the package ctmm (Calabrese et al., 2016) crawl (Johnson et al., 2008),
//' and MScrawl (Michelot and Blackwell, 2019).
//'
//' @name kalman_rcpp
//' @param data Matrix of data, including columns `x`, `y`, `time`, `ID` and `state` (in that order).
//' @param nbStates Integer number of states.
//' @param param Vector of movement parameters (`tau_vel`, `tau_pos`, and `sigma`)
//' @param fixmu Vector of mean locations for the OUF process (`x`, `y`)
//' @param Hmat Matrix of observation error variance (four columns, and one row for each row of data)
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
List kalman_rcpp(arma::mat &data, int nbStates, arma::vec param, arma::vec fixmu, arma::mat &Hmat)
{

  int nbData = data.n_rows;
  int nbID = 0;

  // unpack data
  mat X = data.cols(0, 1);    // x, y
  vec time = data.col(2);     // time
  vec ID = data.col(3);       // ID
  vec S = data.col(4);        // state
  vec dt(nbData, fill::ones); // time intervals
  dt(0) = datum::inf;
  dt.subvec(1, nbData - 1) = diff(time);

  mat aest(4, 3, fill::zeros); // state estimate
  //  x    mu_x    0
  //  vx   mu_vx   0
  //  y    0       mu_y
  //  vy   0       mu_vy
  mat Pest(4, 4, fill::zeros); // covariance estimate

  // unpack parameters
  vec tau_pos = param.subvec(0, nbStates - 1);
  vec tau_vel = param.subvec(nbStates, 2 * nbStates - 1);
  cube sigma(2, 2, nbStates);
  for (int i = 0; i < nbStates; i++)
  {
    if (param.size() / nbStates == 3)
    {
      sigma.slice(i)(0, 0) = param.subvec(2 * nbStates, 3 * nbStates - 1)(i);
      sigma.slice(i)(1, 1) = param.subvec(2 * nbStates, 3 * nbStates - 1)(i);
    }
    else if (param.size() / nbStates == 5)
    {
      sigma.slice(i).diag() = param.subvec(2 * nbStates + (3 * i), 2 * nbStates + (3 * i) + 1);
      sigma.slice(i)(0, 1) = param((2 * nbStates) - 1 + (3 * (i + 1)));
      sigma.slice(i)(1, 0) = param((2 * nbStates) - 1 + (3 * (i + 1)));
    }
  }

  // define all empty matrices and vectors needed for the Kalman Filter and likelihood calculation
  mat Z{{1, 0, 0, 0}, {0, 0, 1, 0}};   // observation model which maps the true state space into the observed space (P, Hk)
  mat I(4, 4, fill::eye);              // identity matrix
  mat H(2, 2, fill::zeros);            // the covariance of the observation noise (error, Rk)
  mat T(4, 4);                         // state transition model which is applied to the previous state xk−1 (Green, Fk)
  mat Q(4, 4);                         // covariance of the process noise (Sigma, Q)
  mat K(4, 2, fill::zeros);            // Kalman Gain ( Gain, Kk )
  cube u(2, 3, nbData, fill::zeros);   // measurement residual (zRes, yk)
  cube iF(2, 2, nbData, fill::zeros);  // inverse of residual covariance, residual "precision" (isRes)
  cube uiF(2, 2, nbData, fill::zeros); // measurement residual * inverse of residual covariance (uisRes)
  cube mu(2, 1, nbData, fill::zeros);
  colvec logdetF(nbData, fill::zeros); // log determinant of residual covariance
  mat mu_out(nbStates, 2);
  mu_out.fill(NA_REAL);

  // vector to keep track of the indices for first positions that are IOU, and size of this vector
  uvec iou(nbData);
  int k = 0;

  // Kalman filter iterations
  for (int i = 0; i < nbData; i++)
  {

    // if first location or new individual
    if (i == 0 || ID(i) != ID(i - 1))
    {
      nbID++;
      // if starting in IOU, keep track of index
      if (std::isinf(tau_pos(S(i) - 1)))
      {
        iou(k) = i;
        k++;
      }
      // reset dt to inf (when new individual)
      dt(i) = datum::inf;

      // initialise state mean
      aest = zeros(4, 3);
      // and initial state covariance matrix
      Pest = makeQ(tau_pos(S(i) - 1), tau_vel(S(i) - 1), sigma.slice(S(i) - 1), dt(i));
    }

    // if starting a new IOU bout, keep track of index
    else if (i > 0 && ID(i) == ID(i - 1) && S(i) != S(i - 1) && std::isinf(tau_pos(S(i) - 1)))
    {
      iou(k) = R_IsNA(X(i, 0)) ? i + 1 : i; // ternary to lookahead of NA inserted by updateState, if doing so
      k++;
    }

    // update our estimate (if, missing obs skip to prediction)
    if (!R_IsNA(X(i, 0)))
    {
      H(0, 0) = Hmat(i, 0);
      H(1, 1) = Hmat(i, 1);
      H(0, 1) = Hmat(i, 2);
      H(1, 0) = Hmat(i, 3);
      mat aobs = join_rows(X.row(i).t(), eye(2, 2));
      // measurement residual (zRes, u)
      u.slice(i) = aobs - (Z * aest);

      // residual covariance (sRes, F)
      mat PestZt = Pest * Z.t();
      PestZt.replace(datum::nan, 0);
      mat ZPestZt = Z * PestZt;
      ZPestZt.replace(datum::nan, 0);
      mat F = ZPestZt + H; // residual covariance (sRes, Sk)
      // iF.slice(i) = F.i();
      iF.slice(i) = PDsolve(F);
      logdetF(i) = det(F) > 0 ? log(std::abs(det(F))) : datum::inf;
      uiF.slice(i) = iF.slice(i) * u.slice(i).submat(0, 1, 1, 2);
      uiF.slice(i).replace(datum::nan, 0);

      // Kalman gain
      K = PestZt * iF.slice(i);
      // if gain is inf or nan, replace with Z
      uvec idx = find_nonfinite(K);
      mat Zt = Z.t();
      K(idx) = Zt(idx);

      // update concurrent state estimate (zCon, aest)
      aest = aest + (K * u.slice(i));
      // update concurrent covariance estimate (sCon, Pest)
      mat J = I - (K * Z);
      mat JPest = J * Pest;
      JPest.replace(datum::nan, 0);
      mat KHKt = K * H * K.t();
      KHKt.replace(datum::nan, 0);
      Pest = JPest * J.t() + KHKt;
    }

    // proceed with forecast
    if (i < nbData - 1 && ID(i) == ID(i + 1))
    {
      T = makeT(tau_pos(S(i + 1) - 1), tau_vel(S(i + 1) - 1), dt(i + 1));
      Q = makeQ(tau_pos(S(i + 1) - 1), tau_vel(S(i + 1) - 1), sigma.slice(S(i + 1) - 1), dt(i + 1));

      // predict state estimate (zFor, aest)
      aest = T * aest;

      // predict covariance estimate (sFor, Pest)
      mat tcon = Pest;
      bool has_inf = tcon.has_inf();
      if (has_inf)
      {
        tcon.replace(datum::inf, 0);
      }

      Pest = T * tcon * T.t() + Q;

      if (has_inf)
      {
        vec Pdiag = Pest.diag();
        bool anyP = Pdiag.has_inf();
        uvec idx(4);
        int l = 0;
        if (anyP)
        {
          for (uword j = 0; j < Pdiag.size(); j++)
          {
            if (std::isinf(Pdiag(j)))
            {
              idx(l) = j;
              l++;
            }
          }
          Pest.rows(idx.head(l)).zeros();
          Pest.cols(idx.head(l)).zeros();
          vec Pdiag2 = Pest.diag();
          Pdiag2(idx.head(l)).fill(datum::inf);
          Pest.diag() = Pdiag2;
        }
      }
    }
  } // end filter iterations

  // calculate state based mu
  // @todo: is element wise multiply on cube slower than reshaping? is there some kind of tensor operation (contraction?)
  for (int i = 0; i < nbStates; i++)
  {
    // if IOU state, use bout start location for mu
    if (std::isinf(tau_pos(i)))
    {
      for (int j = 0; j < nbID; j++)
      {
        uvec iou_idx = intersect(find(S == i + 1), find(ID == j + 1)); // indices for IOU and current individual
        for (uword l = 0; l < iou_idx.size(); l++)
        {
          // find the max iou.head(k) that is not greater than idx
          uword bout_start_idx = find(iou.head(k) <= iou_idx(l)).max();
          mu.slice(iou_idx(l)) = X.row(iou.head(k)(bout_start_idx)).t();
        }
      }
      // otherwise calculate mu from residuals (maximum likelihood solution)
    }
    else
    {
      uvec ouf_idx = intersect(find(S == i + 1), find_finite(X.col(0))); // Find indices where state i occurs
      if (R_IsNA(fixmu(i * 2)) && R_IsNA(fixmu(i * 2 + 1)))
      {
        cube uiF_ouf = uiF.slices(ouf_idx);                           // measurement residual * inverse of residual covariance
        cube u_ouf = u.slices(ouf_idx);                               // measurement residual
        cube D(2, 1, 1);                                              // 'B' in eqns B30, B31
        D.tube(0, 0) = sum(uiF_ouf.tube(0, 0) % u_ouf.tube(0, 0), 2); // data x
        D.tube(1, 0) = sum(uiF_ouf.tube(1, 1) % u_ouf.tube(1, 0), 2); // data y
        cube W = sum(uiF_ouf % u_ouf.tube(0, 1, 1, 2), 2);            // mean x, y
        // mat mu_m = inv(W.slice(0)) * D.slice(0);
        // mat mu_m = solve(W.slice(0), eye(2,2)) * D.slice(0);
        mat mu_m = PDsolve(W.slice(0)) * D.slice(0);
        mu.each_slice(ouf_idx) = mu_m;
        // save mu to output
        mu_out(i, 0) = mu_m(0, 0);
        mu_out(i, 1) = mu_m(1, 0);
      }
      else
      {
        mu.each_slice(ouf_idx) = fixmu.subvec(i * 2, i * 2 + 1);
        // save mu to output
        mu_out(i, 0) = fixmu(i * 2);
        mu_out(i, 1) = fixmu(i * 2 + 1);
      }
    }
  } // end mu

  uvec na_xy = find_nonfinite(X.col(0));

  // detrend & llk
  double llk = 0;
  for (uword i = 0; i < u.n_slices; i++)
  {
    if (!any(i == iou.head(k)) && !any(i == na_xy))
    {
      colvec resid = u.slice(i).col(0) - u.slice(i).cols(1, 2) * mu.slice(i);
      llk -= (logdetF(i) + dot(resid, resid.t() * iF.slice(i))) / 2;
    }
  }

  return List::create(
      Rcpp::Named("llk") = llk,
      Rcpp::Named("mu") = mu_out);
}
