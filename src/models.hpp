#ifndef _MODELS_
#define _MODELS_
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

class Model
{
public:
  virtual Rcpp::NumericMatrix getQ(const double time, const double lat, const double lng) const = 0;
  virtual ~Model() {}
};

std::unique_ptr<Model> createModel(const int &nbStates, const std::string &model_type, const Rcpp::List &params, const double &kappa);

#endif