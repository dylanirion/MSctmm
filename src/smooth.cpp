#include "smooth.hpp"
using namespace Rcpp;
using namespace arma;

//' Kalman smoother
//'
//' This code is adapted from the package ctmm (Calabrese et al., 2016) crawl (Johnson et al., 2008),
//' and MScrawl (Michelot and Blackwell, 2019).
//'
//' @name smooth_rcpp
//' @param data Matrix of data, including columns `x`, `y`, `time`, `ID` and `state` (in that order).
//' @param param Vector of movement parameters (`tau_vel`, `tau_pos`, and `sigma`)
//' @param fixmu Vector of mean locations for the OUF process (`x`, `y`)
//' @param Hmat Matrix of observation error variance (four columns, and one row
//' for each row of data)
//'
//' @return a named List containing the predicted locations and velocities, and variance of these estimates.
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
List smooth_rcpp(const arma::mat &data, const int &nbStates, const arma::vec &param, const arma::vec &fixmu, const arma::mat &Hmat)
{
  int nbData = data.n_rows;

  // unpack data
  mat X = data.cols(0, 1);    // x, y
  vec time = data.col(2);     // time
  vec ID = data.col(3);       // ID
  vec S = data.col(4);        // state
  vec dt(nbData, fill::ones); // time intervals
  dt(0) = datum::inf;
  dt.subvec(1, nbData - 1) = diff(time);

  mat aest(nbData, 4, fill::zeros);     // state estimate
  cube Pest(4, 4, nbData, fill::zeros); // covariance estimate
  mat afor(nbData, 4, fill::zeros);     // state forecast
  cube Pfor(4, 4, nbData, fill::zeros); // covariance forecast

  // unpack parameters
  vec tau_pos = param.subvec(0, nbStates - 1);
  vec tau_vel = param.subvec(nbStates, 2 * nbStates - 1);
  cube sigma(2, 2, nbStates);
  prepare_sigma(param, nbStates, sigma);

  // define all empty matrices and vectors needed for the Kalman Filter and Smoother
  mat Z{{1, 0, 0, 0}, {0, 0, 1, 0}};  // observation model which maps the true state space into the observed space ( P, Hk )
  mat I(4, 4, fill::eye);             // identity matrix
  mat H(2, 2, fill::zeros);           // the covariance of the observation noise ( error, Rk )
  cube T(4, 4, nbData, fill::zeros);  // state transition model which is applied to the previous state xk−1 ( Green, Fk )
  cube Q(4, 4, nbData, fill::zeros);  // covariance of the process noise ( Sigma, Q )
  mat K(4, 2, fill::zeros);           // Kalman Gain ( Gain, Kk )
  mat u(nbData, 2, fill::zeros);      // measurement residual ( zRes, yk )
  cube iF(2, 2, nbData, fill::zeros); // inverse of residual covariance, residual "precision" ( isRes )

  cube L(4, 4, nbData, fill::zeros);

  // forward filter iterations
  for (int i = 0; i < nbData; i++)
  {

    // if first location or new individual
    if (i == 0 || ID(i) != ID(i - 1))
    {
      initialize_state(i, S, tau_pos, tau_vel, sigma, dt, afor.row(i), Pfor.slice(i));
    }

    // update our estimate ( if, missing obs skip to prediction )
    if (!is_observation_missing(X, i))
    {
      H = {{Hmat(i, 0), Hmat(i, 2)}, {Hmat(i, 3), Hmat(i, 1)}};
      // measurement residual (zRes, u)
      u.row(i) = X.row(i) - (Z * afor.row(i).t()).t();

      // residual covariance (sRes, F)
      mat PforZt = Pfor.slice(i) * Z.t();
      PforZt.replace(datum::nan, 0);
      mat ZPforZt = Z * PforZt;
      ZPforZt.replace(datum::nan, 0);
      mat F = ZPforZt + H; // residual covariance ( sRes, Sk )
      // iF.slice(i) = F.i();
      iF.slice(i) = PDsolve(F);

      // Kalman gain
      K = PforZt * iF.slice(i);
      // if gain is inf or nan, replace with Z
      uvec idx = find_nonfinite(K);
      mat Zt = Z.t();
      K(idx) = Zt(idx);

      // update concurrent state estimate (zCon, aest)
      aest.row(i) = afor.row(i) + (K * u.row(i).t()).t();
      // update concurrent covariance estimate (sCon, Pest)
      mat J = I - (K * Z);
      mat JPest = J * Pfor.slice(i);
      JPest.replace(datum::nan, 0);
      mat KHKt = K * H * K.t();
      KHKt.replace(datum::nan, 0);
      Pest.slice(i) = JPest * J.t() + KHKt;
    }

    // proceed with forecast
    if (i < nbData - 1 && ID(i) == ID(i + 1))
    {
      T.slice(i) = makeT(tau_pos(S(i + 1) - 1), tau_vel(S(i + 1) - 1), dt(i + 1));
      Q.slice(i) = makeQ(tau_pos(S(i + 1) - 1), tau_vel(S(i + 1) - 1), sigma.slice(S(i + 1) - 1), dt(i + 1));

      // predict state estimate (zFor, aest)
      afor.row(i + 1) = (T.slice(i) * aest.row(i).t()).t();

      // predict covariance estimate (sFor, Pest)
      mat tcon = Pest.slice(i);
      bool has_inf = tcon.has_inf();
      if (has_inf)
      {
        tcon = tcon.replace(datum::inf, 0);
      }

      Pfor.slice(i + 1) = T.slice(i) * tcon * T.slice(i).t() + Q.slice(i);

      if (has_inf)
      {
        vec Pdiag = Pfor.slice(i + 1).diag();
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
          Pfor.slice(i + 1).rows(idx.head(l)).zeros();
          Pfor.slice(i + 1).cols(idx.head(l)).zeros();
          vec Pdiag2 = Pfor.slice(i + 1).diag();
          Pdiag2(idx.head(l)).fill(datum::inf);
          Pfor.slice(i + 1).diag() = Pdiag2;
        }
      }
    }
  } // end filter iterations

  // backward smoother

  for (uword j = aest.n_rows - 1; j--> 0;)
  {
    if (ID(j) == ID(j + 1))
    {
      mat TL(4, 4);
      if (!is_observation_missing(X, j))
      {
        mat TL = Pest.slice(j) * T.slice(j).t();
      }
      else
      {
        mat TL = Pfor.slice(j) * T.slice(j).t();
      }
      TL.replace(datum::inf, 0);
      // mat INV = solve(Pfor.slice(j + 1), eye(4, 4));
      mat INV = PDsolve(Pfor.slice(j + 1));
      TL = TL * INV;

      vec TLdiag = TL.diag();
      bool has_nan = TLdiag.has_nan();
      uvec idx = find_nonfinite(TLdiag);
      if (has_nan)
      {
        TL.rows(idx).zeros();
        TL.cols(idx).zeros();
        vec TLdiag2 = TL.diag();
        TLdiag2(idx).fill(1);
        TL.diag() = TLdiag2;
      }

      L.slice(j) = TL;
      mat J = I - (L.slice(j) * T.slice(j));
      if (!is_observation_missing(X, j))
      {
        mat JPest = J * Pest.slice(j);
        JPest.replace(datum::nan, 0);
        Pest.slice(j) = JPest * J.t() + (L.slice(j) * (Pest.slice(j + 1) + Q.slice(j)) * L.slice(j).t());
        aest.row(j) = aest.row(j) + (L.slice(j) * (aest.row(j + 1) - afor.row(j + 1)).t()).t();
      }
      else
      {
        // prediction (no est)
        mat JPest = J * Pfor.slice(j);
        JPest.replace(datum::nan, 0);
        Pest.slice(j) = JPest * J.t() + (L.slice(j) * (Pest.slice(j + 1) + Q.slice(j)) * L.slice(j).t());
        aest.row(j) = afor.row(j) + (L.slice(j) * (aest.row(j + 1) - afor.row(j + 1)).t()).t();
      }
    }
    else
    {
      Pest.slice(j) = Pest.slice(j);
      aest.row(j) = aest.row(j);
    }
  }

  return List::create(
      Rcpp::Named("pred") = DataFrame::create(Named("x") = aest.col(0), Named("vx") = aest.col(1), Named("y") = aest.col(2), Named("vy") = aest.col(3)),
      Rcpp::Named("predVar") = Pest);
}
