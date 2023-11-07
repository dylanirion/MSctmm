// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "Q.hpp"
#include "smooth.hpp"

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
 //' Modified from sample_path_mr in ECctmc (Fintzi, 2018)
 // [[Rcpp::export]]
 arma::mat sample_path_mr(const int a, const int b, const double t0, const double t1, const Rcpp::NumericMatrix& Q, const double k) {
   const int limit = 50000;
   const int sublimit = 50000;

   // Get the number of states and initialize vector of states
   int n_states = Q.nrow();
   Rcpp::IntegerVector states = Rcpp::seq_len(n_states);

   // Initialize booleans for whether to keep simulating and whether a
   // valid path has been obtained.
   bool valid_path = false;

   // Initialize objects for storing the sequences of times and states
   std::vector<double> time_vec;
   std::vector<int> state_vec;

   // Sample paths until a valid path has been obtained
   int c = 0;
   int j = 0;
   while(valid_path == false && c < limit) {
     c++;

     // Set boolean to initiate forward sampling
     bool keep_going = true;

     // insert the initial time and state
     time_vec.push_back(t0);
     state_vec.push_back(a);

     // set the current time and state
     Rcpp::NumericVector cur_time(1, t0);
     Rcpp::IntegerVector cur_state(1, a);
     double cur_rate = -Q(cur_state[0] - 1, cur_state[0] - 1);

     // get the state transition probabilities
     Rcpp::NumericVector state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);

     // If the beginning and end states don't match, sample first transition
     if(a != b) {

       // sample the first transition time
       cur_time += -log(1 - Rcpp::runif(1, 0, 1) * (1 - exp(-(t1-t0) * k))) / k;

       // sample the next state
       cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);

       // update the rate of transition out of the new state
       // and update the state transition probabilities
       cur_rate  = -Q(cur_state[0] - 1, cur_state[0] - 1);
       state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);

       // Insert the next state and transition time into the
       // appropriate vectors
       time_vec.push_back(cur_time[0]);
       state_vec.push_back(cur_state[0]);
     }

     // Proceed with forward sampling algorithm
     j = 0;
     while(keep_going == true && j < sublimit) {
       j++;

       if(j == sublimit) {
         keep_going = false;
         valid_path = false;
         time_vec.clear();
         state_vec.clear();
         break;
       }

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
         cur_rate  = -Q(cur_state[0] - 1, cur_state[0] - 1);
         state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);

         // update the state and time vectors
         time_vec.push_back(cur_time[0]);
         state_vec.push_back(cur_state[0]);
       }
     }
   }

   // Add the time and state at the right endpoint
   if(c == limit || (c == 1 && time_vec.size() == 0)) {
     // send auto reject ( fill path with -1s)
     time_vec.clear();
     state_vec.clear();
     time_vec.push_back(-1);
     state_vec.push_back(-1);
   } else{
     // Add the time and state at the right endpoint
     time_vec.push_back(t1);
     state_vec.push_back(b);
   }

   // Collect the sequences of states and times into a matrix
   arma::mat path(time_vec.size(), 2);;
   path.col(0) = arma::conv_to<arma::colvec>::from(time_vec);
   path.col(1) = arma::conv_to<arma::colvec>::from(state_vec);

   return path;
 }


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
arma::mat sample_path_mr2(const int a, const int b, const double t0, const double t1, const double lng0, const double lat0, const double lng1, const double lat1, const int group, const double k, const int nbStates, const arma::vec param, const arma::vec mu, const arma::mat& Hmat, const arma::vec alpha, const arma::vec t_alpha, const String model) {
  const int limit = 1000;

  // initialize vector of states
  Rcpp::IntegerVector states = Rcpp::seq_len(nbStates);

  // Initialize booleans for whether to keep simulating and whether a
  // valid path has been obtained.
  bool valid_path = false;

  // Initialize objects for storing the sequences of times and states
  std::vector<double> time_vec;
  std::vector<int> state_vec;
  std::vector<double> lng_vec;
  std::vector<double> lat_vec;
  std::vector<double> lngvar_vec;
  std::vector<double> latvar_vec;
  std::vector<double> vx_vec;
  std::vector<double> vy_vec;

  // Sample paths until a valid path has been obtained
  int c = 0;
  // TODO: how to memoize getQ?
  while(valid_path == false && c < limit) {
    c++;
    // Set boolean to initiate forward sampling
    bool keep_going = true;

    // insert the initial time and state
    time_vec.push_back(t0);
    state_vec.push_back(a);
    lng_vec.push_back(lng0);
    lat_vec.push_back(lat0);

    // set the current time, state and positions
    Rcpp::NumericVector cur_time(1, t0);
    Rcpp::NumericVector cur_lng(1, lng0);
    Rcpp::NumericVector cur_lat(1, lat0);
    Rcpp::IntegerVector cur_state(1, a);
    double cur_rate;
    Rcpp::NumericMatrix Q;
    try {
      Q = getQ(nbStates, alpha, t_alpha, cur_time[0], cur_lng[0], cur_lat[0], group, model);
      cur_rate = -Q(cur_state[0] - 1, cur_state[0] - 1);
    } catch(std::string e) { // catch sst errors
      if (e.compare("sst") == 0) {
        valid_path = false;
        time_vec.clear();
        state_vec.clear();
        lng_vec.clear();
        lat_vec.clear();
        vx_vec.clear();
        vy_vec.clear();
        break;
      }
    }

    // get the state transition probabilities
    Rcpp::NumericVector state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);

    // If the beginning and end states don't match, sample first transition
    if(a != b) {

      // sample the first transition time
      cur_time += -log(1 - Rcpp::runif(1, 0, 1) * (1 - exp(-(t1-t0) * k))) / k;

      // sample the next state
      cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);

      // predict positions from start and end, with new intermediate transition time
      arma::mat data = {{lng0, lat0, t0, 0, 1.0 * a},
                        {NA_REAL, NA_REAL, cur_time[0], 0, 1.0 * cur_state[0]},
                        {lng1, lat1, t1, 0, 1.0 * b}}; //{x, y, time, ID, state}
      arma::mat newHmat(3, 4);
      newHmat.row(0) = Hmat.row(0);
      newHmat.row(1) = {NA_REAL, NA_REAL, NA_REAL, NA_REAL};
      newHmat.row(2) = Hmat.row(1);
      Rcpp::List smooth = smooth_rcpp(data, nbStates, param, mu, newHmat);
      Rcpp::DataFrame pred = wrap(smooth["pred"]);
      Rcpp::NumericVector x = pred["x"];
      Rcpp::NumericVector y = pred["y"];
      Rcpp::NumericVector vx = pred["vx"];
      Rcpp::NumericVector vy = pred["vy"];
      cur_lng = x[1];
      cur_lat = y[1];


      // update the rate of transition out of the new state
      // and update the state transition probabilities
      try {
        Q = getQ(nbStates, alpha, t_alpha, cur_time[0], cur_lng[0], cur_lat[0], group, model);
        cur_rate  = -Q(cur_state[0] - 1, cur_state[0] - 1);
      } catch(std::string e) { // catch sst errors
        if (e.compare("sst") == 0) {
          valid_path = false;
          time_vec.clear();
          state_vec.clear();
          lng_vec.clear();
          lat_vec.clear();
          vx_vec.clear();
          vy_vec.clear();
          break;
        }
      }
      state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);

      // Insert the next state and transition time into the
      // appropriate vectors
      time_vec.push_back(cur_time[0]);
      state_vec.push_back(cur_state[0]);
      lng_vec.push_back(cur_lng[0]);
      lat_vec.push_back(cur_lat[0]);
      vx_vec.push_back(vx[1]);
      vy_vec.push_back(vy[1]);
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
          lng_vec.clear();
          lat_vec.clear();
          vx_vec.clear();
          vy_vec.clear();
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
          lng_vec.clear();
          lat_vec.clear();
          vx_vec.clear();
          vy_vec.clear();
        }

        // If the next transition occurs before the right
        // endpoint, sample the next state
      } else {
        // sample the next state
        cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);

        // predict positions from start and end, with new intermediate transition time
        arma::mat data(time_vec.size() + 2, 5);
        arma::mat newHmat(time_vec.size() + 2, 4);
        data.row(0) = {lng0, lat0, t0, 0, 1.0 * a};
        newHmat.row(0) = Hmat.row(0);
        for(unsigned i = 1; i < time_vec.size(); i++) {
          data.row(i) = { lng_vec[i], lat_vec[i], time_vec[i], 0, 1.0 * state_vec[i]};
          newHmat.row(i) = {vx_vec[i], vy_vec[i], 0, 0};
        }
        data.row(time_vec.size()) = {NA_REAL, NA_REAL, cur_time[0], 0, 1.0 * cur_state[0]};
        newHmat.row(time_vec.size()) = {NA_REAL, NA_REAL, NA_REAL, NA_REAL};
        newHmat.row(time_vec.size() + 1) = Hmat.row(1);
        data.row(time_vec.size() + 1) = {lng1, lat1, t1, 0, 1.0 * b};
        Rcpp::List smooth = smooth_rcpp(data, nbStates, param, mu, newHmat);
        Rcpp::DataFrame pred = wrap(smooth["pred"]);
        Rcpp::NumericVector x = pred["x"];
        Rcpp::NumericVector y = pred["y"];
        Rcpp::NumericVector vx = pred["vx"];
        Rcpp::NumericVector vy = pred["vy"];
        cur_lng = x[time_vec.size() + 1];
        cur_lat = y[time_vec.size() + 1];

        // update the rate of transition out of the new state
        // and update the state transition probabilities
        try {
          Q = getQ(nbStates, alpha, t_alpha, cur_time[0], cur_lng[0], cur_lat[0], group, model);
          cur_rate  = -Q(cur_state[0] - 1, cur_state[0] - 1);
        } catch(std::string e) { // catch sst errors
          if (e.compare("sst") == 0) {
            valid_path = false;
            time_vec.clear();
            state_vec.clear();
            lng_vec.clear();
            lat_vec.clear();
            vx_vec.clear();
            vy_vec.clear();
            break;
          }
        }
        state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);
        // update the state and time vectors
        time_vec.push_back(cur_time[0]);
        state_vec.push_back(cur_state[0]);
        lng_vec.push_back(cur_lng[0]);
        lat_vec.push_back(cur_lat[0]);
        vx_vec.push_back(vx[time_vec.size() + 1]);
        vy_vec.push_back(vy[time_vec.size() + 1]);
      }
    }
  }

  if(c == limit || (c == 1 && time_vec.size() == 0)) {
    // send auto reject ( fill path with -1s)
    time_vec.push_back(-1);
    state_vec.push_back(-1);
  } else{
    // Add the time and state at the right endpoint
    time_vec.push_back(t1);
    state_vec.push_back(b);
  }

  // Collect the sequences of states and times into a matrix
  arma::mat path(time_vec.size(), 2);
  path.col(0) = arma::conv_to<arma::colvec>::from(time_vec);
  path.col(1) = arma::conv_to<arma::colvec>::from(state_vec);

  return path;
}
