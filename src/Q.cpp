#include <RcppArmadillo.h>
#include <ctime>
#include <iomanip>
#include <proj.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

PJ_CONTEXT *C = proj_context_create();
PJ *P = proj_create_crs_to_crs(C, "+proj=tpeqd +lat_1=-34.09303 +lon_1=22.19052 +lat_2=-34.09303 +lon_2=22.25142 +x_0=0 +y_0=0 +datum=WGS84", "EPSG:4326", NULL);

// Make Q matrix
// [[Rcpp::export]]
Rcpp::NumericMatrix getQ(const int nbStates, arma::vec alpha, arma::vec t_alpha, const time_t time, const double lng, const double lat, const String model) {
  tm *t = localtime(&time);
  int yday = t->tm_yday;

  PJ_COORD input_coords, output_coords; // https://proj.org/development/reference/datatypes.html#c.PJ_COORD
  input_coords = proj_coord(lat, lng, 0, 0);
  output_coords = proj_trans(P, PJ_FWD, input_coords);

  mat Q(nbStates, nbStates, fill::ones);
  Q.diag() = Q.diag() * -1;

  if (model == "time_out_time_in") {
    // time-varying rate in and out
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
  } else if (model == "sst_out_sst_in") {
    //TODO HOW TO HANDLE POINTS ON LAND OR MISSING SST?
    Rcpp::Environment rerddap = Rcpp::Environment::namespace_env("rerddap");
    Rcpp::Function f = rerddap["griddap"];
    Rcpp::CharacterVector times = Rcpp::CharacterVector::create(std::to_string(t->tm_year + 1900) + "-" + std::to_string(t->tm_mon + 1) + "-" + std::to_string(t->tm_mday), std::to_string(t->tm_year + 1900) + "-" + std::to_string(t->tm_mon + 1) + "-" + std::to_string(t->tm_mday));
    Rcpp::NumericVector lngs = Rcpp::NumericVector::create(output_coords.xy.x, output_coords.xy.x);
    Rcpp::NumericVector lats = Rcpp::NumericVector::create(output_coords.xy.y, output_coords.xy.y);
    Rcpp::DataFrame sst_df = f("jplMURSST41", Rcpp::Named("time", times), Rcpp::Named("longitude", lngs), Rcpp::Named("latitude", lats), Rcpp::Named("fields", "analysed_sst"), Rcpp::Named("fmt", "csv"));
    Rcpp::NumericVector sst = sst_df["analysed_sst"];
    Rcout << std::to_string(t->tm_year + 1900) + "-" + std::to_string(t->tm_mon + 1) + "-" + std::to_string(t->tm_mday) + " " + std::to_string(t->tm_hour) + ":" + std::to_string(t->tm_min) + ":" + std::to_string(t->tm_sec) << ", " << output_coords.xy.x << ", " << output_coords.xy.y << ", " << sst[0] << endl;
    // sst-varying rate in and out
    //FB -> trans
    Q(0,4) = alpha(0)/(1+exp(-alpha(0) * (sst[0] - t_alpha(0))));
    Q(0,0) = Q(0,4) * -1;
    //GB -> trans
    Q(1,4) = alpha(0)/(1+exp(-alpha(0) * (sst[0] - t_alpha(0))));
    Q(1,1) = Q(1,4) * -1;
    //MB -> trans
    Q(2,4) = alpha(0)/(1+exp(-alpha(0) * (sst[0] - t_alpha(0))));
    Q(2,2) = Q(2,4) * -1;
    //AB -> trans
    Q(3,4) = alpha(0)/(1+exp(-alpha(0) * (sst[0] - t_alpha(0))));
    Q(3,3) = Q(3,4) * -1;
    // trans -> FB
    Q(4,0) = alpha(1)/(1+exp(-alpha(1) * (sst[0] - t_alpha(1))))/4;
    // trans -> GB
    Q(4,1) = alpha(1)/(1+exp(-alpha(1) * (sst[0] - t_alpha(1))))/4;
    // trans -> MB
    Q(4,2) = alpha(1)/(1+exp(-alpha(1) * (sst[0] - t_alpha(1))))/4;
    // trans -> AB
    Q(4,3) = alpha(1)/(1+exp(-alpha(1) * (sst[0] - t_alpha(1))))/4;
    Q(4,4) = -(alpha(1)/(1+exp(-alpha(1) * (sst[0] - t_alpha(1)))));
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
