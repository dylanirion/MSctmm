// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
// #include "Q.hpp"
#include "models.hpp"
#include "smooth.hpp"
using namespace arma;
using namespace Rcpp;

mat create_path(const std::vector<double> &time_vec, const std::vector<int> &state_vec);

mat create_invalid_path();

mat create_data_matrix(const std::vector<double> &time_vec, const std::vector<int> &state_vec,
                             const std::vector<double> &lng_vec, const std::vector<double> &lat_vec,
                             const double &t1, const int &a, const int &b, const double &cur_time, const int &pot_state);

mat create_hmat(const arma::mat &Hmat, const size_t &size);

//' Simulate a sample path from an endpoint conditioned CTMC by modified
//' rejection sampling.
//'
//' @param a,b States at the interval endpoints, provided as integers
//'    corresponding to rows of the CTMC rate matrix.
//' @param t0,t1 times of the interval endpoints
//'
//' @return matrix whose first column is the sequence of transition times
//' bookended by interval endpoints, and whose second column is the sequence of
//' states
//'
//' Modified from sample_path_mr in ECctmc (Fintzi, 2018)
// [[Rcpp::export]]
arma::mat sample_path_mr_(const int a, const int b, const double t0, const double t1, const double lng0, const double lat0, const double lng1, const double lat1, const int group, const double k, const int &nbStates, const arma::vec &param, const arma::vec &mu, const arma::mat &Hmat, Rcpp::List rateparam, const String model)
{
  // initialize vector of states
  IntegerVector states = seq_len(nbStates);

  // initialize objects for storing the sequences of times and states
  std::vector<double> time_vec, lng_vec, lat_vec, vx_vec, vy_vec;
  std::vector<int> state_vec;

  //if non-na mus, choose spatial model, (and add mu to rateparam)
  // was previously only doing this for SpatialNA model but need it elsewhere and it doesn't hurt
  uvec finite = arma::find_finite(mu);
  //if (model == "NA" && finite.size() > 1) {
  if (finite.size() > 1) {
    rateparam.push_back(wrap(mu), "mu");
  }
  auto mod = createModel(nbStates, (model == "NA" && finite.size() > 1) ? "SpatialNA" : model, rateparam, k);

  // sample paths until a valid path has been obtained
  for (int i = 0; i < 100; i++)
  {
    // insert the initial time and state
    time_vec = {t0};
    state_vec = {a};
    lng_vec = {lng0};
    lat_vec = {lat0};
    vx_vec.clear();
    vy_vec.clear();

    // set the current time, state and positions
    double cur_time = t0;
    int cur_state = a;
    double cur_lng = lng0;
    double cur_lat = lat0;

    // ECctmc had something like this, but I feel this incorrectly precludes
    // switching over long time gaps,
    // or incorrectly specified init state sequences?
    // e.g. a sequence entirely of state "a" will never switch
    //if (a == b)
    //{
    //  time_vec.push_back(t1);
    //  state_vec.push_back(b);
    //  return create_path(time_vec, state_vec);
    //}

    NumericMatrix Q;
    try
    {
      Q = mod->getQ(cur_time, cur_lng, cur_lat);
    }
    catch (const std::string &e)
    { // catch sst errors
      if (e == "sst")
      {
        break;
      }
    }
    double cur_rate = -Q(cur_state - 1, cur_state - 1);

    // get the state transition probabilities
    NumericVector state_probs = pmax(Q(cur_state - 1, _), 0);

    int j = 0;
    // Proceed with forward sampling algorithm
    while (cur_time < t1 && j < 2000)
    {
      j++;
      // check if the state is an absorbing state
      if (is_true(all(state_probs == 0)))
      {
        if (cur_state == b)
        {
          time_vec.push_back(t1);
          state_vec.push_back(b);
          return create_path(time_vec, state_vec);
        }
        break;
      }

      if (cur_time == t0)
      {
        // sample the first potential transition time
        cur_time += -log(1 - runif(1, 0, 1)[0] * (1 - exp(-(t1 - t0) * k))) / k;
      }
      else
      {
        // Sample the next potential transition time
        cur_time += rexp(1, cur_rate)[0];
      }

      // If the next time is after the right endpoint, stop sampling
      if (cur_time > t1)
      {
        // and determine if the path is valid
        if (cur_state == b)
        {
          time_vec.push_back(t1);
          state_vec.push_back(b);
          return create_path(time_vec, state_vec);
        }
        break;
      }

      // sample the next potential state
      int pot_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs)[0];

      // predict positions from start and end, with new intermediate transition time
      mat data = create_data_matrix(time_vec, state_vec, lng_vec, lat_vec, t1, a, b, cur_time, pot_state);
      mat newHmat = create_hmat(Hmat, time_vec.size());
      List smooth = smooth_rcpp(data, nbStates, param, mu, newHmat);
      DataFrame pred = wrap(smooth["pred"]); //TODO: REDUNDANT WRAPPING

      double pot_lng = as<NumericVector>(pred["x"])[time_vec.size()]; //TODO: REDUNDANT CASTING
      double pot_lat = as<NumericVector>(pred["y"])[time_vec.size()]; //TODO: REDUNDANT CASTING

      // update the rate of transition out of the new state
      // and update the state transition probabilities
      try
      {
        Q = mod->getQ(cur_time, pot_lng, pot_lat);
        double pot_rate = -Q(pot_state - 1, pot_state - 1);

        // if actual switch, update time, state, rate and transition probabilities
        if (runif(1)[0] < (pot_rate / k))
        {
          state_probs = pmax(Q(cur_state - 1, _), 0);
          cur_state = pot_state;
          cur_rate = pot_rate;
          cur_lng = pot_lng;
          cur_lat = pot_lat;

          // Insert the next state and transition time into the
          // appropriate vectors
          time_vec.push_back(cur_time);
          state_vec.push_back(cur_state);
          lng_vec.push_back(cur_lng);
          lat_vec.push_back(cur_lat);
          vx_vec.push_back(as<NumericVector>(pred["vx"])[time_vec.size()]); //TODO: REDUNDANT CASTING
          vy_vec.push_back(as<NumericVector>(pred["vy"])[time_vec.size()]); //TODO: REDUNDANT CASTING
        }
      }
      catch (const std::string &e)
      {
        if (e == "sst")
        {
          break;
        }
      }
    }
  }

  // If we've reached here, we've failed to generate a valid path
  // send auto reject (fill path with -1s)
  return create_invalid_path();
}

mat create_path(const std::vector<double> &time_vec, const std::vector<int> &state_vec)
{
  // Collect the sequences of states and times into a matrix
  mat path(time_vec.size(), 2);
  path.col(0) = conv_to<colvec>::from(time_vec);
  path.col(1) = conv_to<colvec>::from(state_vec);
  return path;
}

mat create_invalid_path()
{
  mat path = {{-1, -1}};
  return path;
}

mat create_data_matrix(const std::vector<double> &time_vec, const std::vector<int> &state_vec,
                             const std::vector<double> &lng_vec, const std::vector<double> &lat_vec,
                             const double &t1, const int &a, const int &b, const double &cur_time, const int &pot_state)
{
  mat data(time_vec.size() + 2, 5);
  for (size_t i = 0; i < time_vec.size(); ++i)
  {
    data.row(i) = {lng_vec[i], lat_vec[i], time_vec[i], 0, static_cast<double>(state_vec[i])};
  }
  data.row(time_vec.size()) = {NA_REAL, NA_REAL, cur_time, 0, static_cast<double>(pot_state)};
  data.row(time_vec.size() + 1) = {lng_vec.back(), lat_vec.back(), t1, 0, static_cast<double>(b)};
  return data;
}

mat create_hmat(const mat &Hmat, const size_t &size)
{
  mat newHmat(size + 2, 4, fill::none);
  newHmat.row(0) = Hmat.row(0);
  newHmat.rows(1, size) = mat(size, 4, fill::none);
  newHmat.row(size + 1) = Hmat.row(1);
  return newHmat;
}
