// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "Q.hpp"

using namespace arma;
using namespace Rcpp;

//' Simulate a sample path from an endpoint conditioned CTMC by modified
//' rejection sampling.
//'
//' @param a,b States at the interval endpoints, provided as integers
//'    corresponding to rows of the CTMC rate matrix.
//' @param t0,t1 times of the interval endpoints
//' @param Q CTMC rate matrix
//'
//' @return matrix whose first column is the sequence of transition times
//' bookended by interval endpoints, and whose second column is the sequence of
//' states
//'
//' Modified from sample_math_mr in ECctmc (Fintzi, 2018)
// [[Rcpp::export]]
arma::mat sample_path_mr2(const int a, const int b, const double t0, const double t1, const double k, const int nbStates, arma::vec alpha, arma::vec t_alpha, String model) {

  // initialize vector of states
  Rcpp::IntegerVector states = Rcpp::seq_len(nbStates);

  // Initialize booleans for whether to keep simulating and whether a
  // valid path has been obtained.
  bool valid_path = false;

  // Initialize objects for storing the sequences of times and states
  std::vector<double> time_vec;
  std::vector<int> state_vec;

  // Sample paths until a valid path has been obtained
  while(valid_path == false) {

    // Set boolean to initiate forward sampling
    bool keep_going = true;

    // insert the initial time and state
    time_vec.push_back(t0);
    state_vec.push_back(a);
    double lng = 0;
    double lat = 0;

    // set the current time and state
    Rcpp::NumericVector cur_time(1, t0);
    Rcpp::IntegerVector cur_state(1, a);
    double cur_rate = -getQ(nbStates, alpha, t_alpha, cur_time[0], lng, lat, model)(cur_state[0] - 1, cur_state[0] - 1);

    // get the state transition probabilities
    Rcpp::NumericVector state_probs = pmax(getQ(nbStates, alpha, t_alpha, cur_time[0], lng, lat, model)(cur_state[0] - 1, _ ), 0);

    // If the beginning and end states don't match, sample first transition
    if(a != b) {

      // sample the first transition time
      cur_time += -log(1 - Rcpp::runif(1, 0, 1) * (1 - exp(-(t1-t0) * k))) / k;

      // sample the next state
      cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);

      // update the rate of transition out of the new state
      // and update the state transition probabilities
      cur_rate  = -getQ(nbStates, alpha, t_alpha, cur_time[0], lng, lat, model)(cur_state[0] - 1, cur_state[0] - 1);
      state_probs = pmax(getQ(nbStates, alpha, t_alpha, cur_time[0], lng, lat, model)(cur_state[0] - 1, _ ), 0);

      // Insert the next state and transition time into the
      // appropriate vectors
      time_vec.push_back(cur_time[0]);
      state_vec.push_back(cur_state[0]);
    }

    // Proceed with forward sampling algorithm
    while(keep_going == true) {

      // check if the state is an absorbing state
      if(is_true(all(state_probs == 0))) {

        // stop sampling forward
        keep_going = false;

        if(cur_state[0] == b) {
          valid_path = true;
        } else {
          valid_path = false;
          time_vec.clear();
          state_vec.clear();
        }

        break;
      }

      // Sample the next transition time
      cur_time += Rcpp::rexp(1, cur_rate);

      // If the next time is after the right endpoint, stop
      // sampling and determine if the path is valid
      if(cur_time[0] > t1) {

        // Stop forward sampling
        keep_going = false;

        // Determine if the path is valid
        if(cur_state[0] == b) {
          valid_path = true;

        } else {
          valid_path = false;
          time_vec.clear();
          state_vec.clear();
        }

        // If the next transition occurs before the right
        // endpoint, sample the next state
      } else {
        // sample the next state
        cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);

        // update the rate of transition out of the new state
        // and update the state transition probabilities
        cur_rate  = -getQ(nbStates, alpha, t_alpha, cur_time[0], lng, lat, model)(cur_state[0] - 1, cur_state[0] - 1);
        state_probs = pmax(getQ(nbStates, alpha, t_alpha, cur_time[0], lng, lat, model)(cur_state[0] - 1, _ ), 0);
        // update the state and time vectors
        time_vec.push_back(cur_time[0]);
        state_vec.push_back(cur_state[0]);
      }
    }
  }

  // Add the time and state at the right endpoint
  time_vec.push_back(t1);
  state_vec.push_back(b);

  // Collect the sequences of states and times into a matrix
  arma::mat path(time_vec.size(), 2);;
  path.col(0) = arma::conv_to<arma::colvec>::from(time_vec);
  path.col(1) = arma::conv_to<arma::colvec>::from(state_vec);

  return path;
}
