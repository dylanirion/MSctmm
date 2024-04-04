#include <RcppArmadillo.h>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <proj.h>
#include <chrono>
#include <thread>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

PJ_CONTEXT *C = proj_context_create();
PJ *P = proj_create_crs_to_crs(C, "+proj=tpeqd +lat_1=-34.09303 +lon_1=22.19052 +lat_2=-34.09303 +lon_2=22.25142 +x_0=0 +y_0=0 +datum=WGS84", "EPSG:4326", NULL);

double clamp(double d, double min, double max) {
  const double t = d < min ? min : d;
  return t > max ? max : t;
}

Rcpp::Environment rerddap = Rcpp::Environment::namespace_env("rerddap");
//Rcpp::Function c = rerddap["cache_setup"];
//SEXP tmp = c(Rcpp::Named("full_path", "/home/dylan/rerddap/"));
Rcpp::Function f = rerddap["griddap"];


// Make Q matrix
// [[Rcpp::export]]
Rcpp::NumericMatrix getQ(const int nbStates, arma::vec alpha, arma::vec t_alpha, const time_t time, const double lng, const double lat, const int group, const String model) {
  struct tm t = *localtime(&time);
  int yday = t.tm_yday;

  PJ_COORD input_coords, output_coords; // https://proj.org/development/reference/datatypes.html#c.PJ_COORD
  input_coords = proj_coord(lat, lng, 0, 0);
  output_coords = proj_trans(P, PJ_FWD, input_coords);

  mat Q(nbStates, nbStates, fill::ones);
  Q.diag() = Q.diag() * -1;

  if (model == "time_out_time_in") {
    // time-varying rate in and out
    t_alpha(0) = clamp(t_alpha(0), 0, 365);
    t_alpha(1) = clamp(t_alpha(1), 0, 365);
    //FB -> trans
    Q(0,4) = alpha(0)/(1+exp(-alpha(0) * (yday - t_alpha(0))));
    Q(0,0) = Q(0,4) * -1;
    //GB -> trans
    Q(1,4) = alpha(0)/(1+exp(-alpha(0) * (yday - t_alpha(0))));
    Q(1,1) = Q(1,4) * -1;
    //MB -> trans
    Q(2,4) = alpha(0)/(1+exp(-alpha(0) * (yday - t_alpha(0))));
    Q(2,2) = Q(2,4) * -1;
    //AB -> trans
    Q(3,4) = alpha(0)/(1+exp(-alpha(0) * (yday - t_alpha(0))));
    Q(3,3) = Q(3,4) * -1;
    // trans -> FB
    Q(4,0) = alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1))))/4;
    // trans -> GB
    Q(4,1) = alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1))))/4;
    // trans -> MB
    Q(4,2) = alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1))))/4;
    // trans -> AB
    Q(4,3) = alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1))))/4;
    Q(4,4) = -(alpha(1)/(1+exp(-alpha(1) * (yday - t_alpha(1)))));
    //Impossible transitions
    Q(0,1) = 0; // This might actually be possible (GB->FB)
    Q(0,2) = 0;
    Q(0,3) = 0;
    Q(1,0) = 0; // This might actually be possible (FB->GB)
    Q(1,2) = 0;
    Q(1,3) = 0;
    Q(2,0) = 0;
    Q(2,1) = 0;
    Q(2,3) = 0;
    Q(3,0) = 0;
    Q(3,1) = 0;
    Q(3,2) = 0;
  } else if (model == "time_out_time_in_group") {
    // time-varying rate in and out, with n group-specific rates (first n rates are out, next n rates are in)
    // (this actually functions identically to above)
    for(unsigned i = 0; i < t_alpha.size(); i++) {
      t_alpha(i) = clamp(t_alpha(i), 0, 365);
    }
    //FB -> trans
    Q(0,4) = alpha(group - 1)/(1+exp(-alpha(group - 1) * (yday - t_alpha(group - 1))));
    Q(0,0) = Q(0,4) * -1;
    //GB -> trans
    Q(1,4) = alpha(group - 1)/(1+exp(-alpha(group - 1) * (yday - t_alpha(group - 1))));
    Q(1,1) = Q(1,4) * -1;
    //MB -> trans
    Q(2,4) = alpha(group - 1)/(1+exp(-alpha(group - 1) * (yday - t_alpha(group - 1))));
    Q(2,2) = Q(2,4) * -1;
    //AB -> trans
    Q(3,4) = alpha(group - 1)/(1+exp(-alpha(group - 1) * (yday - t_alpha(group - 1))));
    Q(3,3) = Q(3,4) * -1;
    // trans -> FB
    Q(4,0) = alpha((t_alpha.size() / 2) + (group - 1))/(1+exp(-alpha((t_alpha.size() / 2) + (group - 1)) * (yday - t_alpha((t_alpha.size() / 2) + (group - 1)))))/4;
    // trans -> GB
    Q(4,1) = alpha((t_alpha.size() / 2) + (group - 1))/(1+exp(-alpha((t_alpha.size() / 2) + (group - 1)) * (yday - t_alpha((t_alpha.size() / 2) + (group - 1)))))/4;
    // trans -> MB
    Q(4,2) = alpha((t_alpha.size() / 2) + (group - 1))/(1+exp(-alpha((t_alpha.size() / 2) + (group - 1)) * (yday - t_alpha((t_alpha.size() / 2) + (group - 1)))))/4;
    // trans -> AB
    Q(4,3) = alpha((t_alpha.size() / 2) + (group - 1))/(1+exp(-alpha((t_alpha.size() / 2) + (group - 1)) * (yday - t_alpha((t_alpha.size() / 2) + (group - 1)))))/4;
    Q(4,4) = -(alpha((t_alpha.size() / 2) + (group - 1))/(1+exp(-alpha((t_alpha.size() / 2) + (group - 1)) * (yday - t_alpha((t_alpha.size() / 2) + (group - 1))))));
    //Impossible transitions
    Q(0,1) = 0; // This might actually be possible (GB->FB)
    Q(0,2) = 0;
    Q(0,3) = 0;
    Q(1,0) = 0; // This might actually be possible (FB->GB)
    Q(1,2) = 0;
    Q(1,3) = 0;
    Q(2,0) = 0;
    Q(2,1) = 0;
    Q(2,3) = 0;
    Q(3,0) = 0;
    Q(3,1) = 0;
    Q(3,2) = 0;
  } else if (model == "sst_out_sst_in") {
    //Checking if point falls on land will make this already slow lookup even longer
    //hacky solution is to expand a point outwards by 0.1 degrees until we get sst values and take the average
    //I do this for a maximum of 10 iterations before we bail out and consider it a bad simulation
    time_t curr_time;
    tm * curr_tm;
    char time_string[100];
    std::time(&curr_time);
    curr_tm = localtime(&curr_time);
    strftime(time_string, 50, "%T", curr_tm);
    int i = 0;
    bool not_valid = true;
    double sst = 0.0;
    while(not_valid) {
      Rcout << "\33[2K\r" << time_string << " Q" << "(" << i << ") " << (round(output_coords.xy.x * 100) / 100)<< ", " << (round(output_coords.xy.y * 100) / 100) << ": " << (not_valid ? "not valid" : "valid");
      std::ostringstream oss;
      oss << (t.tm_year + 1900) << "-" << std::setw(2) << std::setfill('0') << (t.tm_mon + 1) << "-" << std::setw(2) << std::setfill('0') << t.tm_mday;
      std::string date = oss.str();
      Rcpp::CharacterVector times = Rcpp::CharacterVector::create(date, date);
      Rcpp::NumericVector lngs = Rcpp::NumericVector::create((round(output_coords.xy.x * 100) / 100) - (i * 0.1), (round(output_coords.xy.x * 100) / 100) + (i * 0.1));
      Rcpp::NumericVector lats = Rcpp::NumericVector::create((round(output_coords.xy.y * 100) / 100) - (i * 0.1), (round(output_coords.xy.y * 100) / 100) + (i * 0.1));
      Rcpp::DataFrame sst_df;
      try {
        sst_df = f("jplMURSST41", Rcpp::Named("time", times), Rcpp::Named("longitude", lngs), Rcpp::Named("latitude", lats), Rcpp::Named("fields", "analysed_sst"), Rcpp::Named("fmt", "csv"), Rcpp::Named("url", "https://coastwatch.pfeg.noaa.gov/erddap"), Rcpp::Named("callopts", Rcpp::List::create(Rcpp::Named("ipresolve") = 1)));
      } catch(...) {
        Rcout << "\33[2K\r" << time_string << " Q" << "(" << i << ") " << (round(output_coords.xy.x * 100) / 100) << ", " << (round(output_coords.xy.y * 100) / 100) << ": " << (not_valid ? "not valid" : "valid") << " (waiting)";
        //std::this_thread::sleep_for(std::chrono::minutes(1));
        continue;
      }
      Rcpp::NumericVector sst_vec = sst_df["analysed_sst"];
      not_valid = sum(!is_nan(sst_vec)) == 0;
      if (not_valid) {
        if(i == 10) {
          throw std::string("sst");
        }
        i++;
        //std::this_thread::sleep_for(std::chrono::milliseconds(10 * i)); //sleep for 0.01 seconds, increasing as we try more and more
      } else {
        Rcpp::NumericVector temp = sst_vec[!is_nan(sst_vec)];
        sst = mean(temp);
      }
    }
    Rcout << "\33[2K\r" << time_string << " Q" << "(" << i << ") " << (round(output_coords.xy.x * 100) / 100) << ", " << (round(output_coords.xy.y * 100) / 100) << ": " << (not_valid ? "not valid" : "valid");
    // sst-varying rate in and out
    //FB -> trans
    Q(0,4) = alpha(0)/(1+exp(-alpha(0) * (sst - t_alpha(0))));
    Q(0,0) = Q(0,4) * -1;
    //GB -> trans
    Q(1,4) = alpha(0)/(1+exp(-alpha(0) * (sst - t_alpha(0))));
    Q(1,1) = Q(1,4) * -1;
    //MB -> trans
    Q(2,4) = alpha(0)/(1+exp(-alpha(0) * (sst - t_alpha(0))));
    Q(2,2) = Q(2,4) * -1;
    //AB -> trans
    Q(3,4) = alpha(0)/(1+exp(-alpha(0) * (sst - t_alpha(0))));
    Q(3,3) = Q(3,4) * -1;
    // trans -> FB
    Q(4,0) = alpha(1)/(1+exp(-alpha(1) * (sst - t_alpha(1))))/4;
    // trans -> GB
    Q(4,1) = alpha(1)/(1+exp(-alpha(1) * (sst - t_alpha(1))))/4;
    // trans -> MB
    Q(4,2) = alpha(1)/(1+exp(-alpha(1) * (sst - t_alpha(1))))/4;
    // trans -> AB
    Q(4,3) = alpha(1)/(1+exp(-alpha(1) * (sst - t_alpha(1))))/4;
    Q(4,4) = -(alpha(1)/(1+exp(-alpha(1) * (sst - t_alpha(1)))));
    //Impossible transitions
    Q(0,1) = 0; // This might actually be possible (GB->FB)
    Q(0,2) = 0;
    Q(0,3) = 0;
    Q(1,0) = 0; // This might actually be possible (FB->GB)
    Q(1,2) = 0;
    Q(1,3) = 0;
    Q(2,0) = 0;
    Q(2,1) = 0;
    Q(2,3) = 0;
    Q(3,0) = 0;
    Q(3,1) = 0;
    Q(3,2) = 0;
  }

  return Rcpp::wrap(Q);
}
