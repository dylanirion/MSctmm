#include "models.hpp"

// TODO: could add a needs_position flag to Model and check in simulation to prevent unnecessary smoothing
// TODO: group and other parameters necessary for models?

// Null Model
// In this transition rate model, all animals share the same rate generator matrix
// defined by rates out of each resident state, and rates into each resident state
// Some transitions are impossible
class NullModel : public Model
{
private:
  int nbStates;
  // vector of rates, first nbStates - 1 are rates out of resident states, second nbStates - 1 are rates into resident states
  Rcpp::NumericVector alpha;
  arma::mat Q;

public:
  NullModel(const int nbStates, const Rcpp::NumericVector alpha)
      : nbStates(nbStates),
        alpha(alpha),
        Q(nbStates, nbStates, arma::fill::ones)
  {
    //TODO: do this with a loop and automatically detect range states 
    Q.diag() = Q.diag() * -1;
    // FB -> trans
    Q(0, 4) = alpha(0);
    Q(0, 0) = Q(0, 4) * -1;
    // GB -> trans
    Q(1, 4) = alpha(1);
    Q(1, 1) = Q(1, 4) * -1;
    // MB -> trans
    Q(2, 4) = alpha(2);
    Q(2, 2) = Q(2, 4) * -1;
    // AB -> trans
    Q(3, 4) = alpha(3);
    Q(3, 3) = Q(3, 4) * -1;
    // trans -> FB
    Q(4, 0) = alpha(4);
    // trans -> GB
    Q(4, 1) = alpha(5);
    // trans -> MB
    Q(4, 2) = alpha(6);
    // trans -> AB
    Q(4, 3) = alpha(7);
    Q(4, 4) = (Q(4, 0) + Q(4, 1) + Q(4, 2) + Q(4, 3)) * -1;
    // Set impossible transitions
    Q(0, 1) = 0; // This might actually be possible (GB->FB)
    Q(0, 2) = 0;
    Q(0, 3) = 0;
    Q(1, 0) = 0; // This might actually be possible (FB->GB)
    Q(1, 2) = 0;
    Q(1, 3) = 0;
    Q(2, 0) = 0;
    Q(2, 1) = 0;
    Q(2, 3) = 0;
    Q(3, 0) = 0;
    Q(3, 1) = 0;
    Q(3, 2) = 0;
  }

  Rcpp::NumericMatrix getQ(const double time, const double lat, const double lng) const override
  {
    return Rcpp::wrap(Q);
  }
};

// NA Model
// In this transition rate model, individuals have their own rate generator matrices
// derived from Gibbs sampling on the number of intervals, the time spent in each state
// and the number of transitions from their current state sequence
class NAModel : public Model
{
private:
  int nbStates;
  //NB times and states of switches only! not full data
  std::vector<double> times;
  std::vector<int> states;
  Rcpp::NumericVector priorShape;
  Rcpp::NumericVector priorRate;
  Rcpp::NumericVector priorCon;
  double kappa;
  arma::mat Q;

public:
  NAModel(const int nbStates, const Rcpp::DataFrame data, const Rcpp::NumericVector priorShape, const Rcpp::NumericVector priorRate, const Rcpp::NumericVector priorCon, const double kappa)
      : nbStates(nbStates),
        times(Rcpp::as<std::vector<double>>(data["time"])),
        states(Rcpp::as<std::vector<int>>(data["state"])),
        priorShape(priorShape),
        priorRate(priorRate),
        priorCon(priorCon),
        kappa(kappa),
        Q(nbStates, nbStates, arma::fill::zeros)
  {
    std::vector<double> timeInStates(nbStates, 0.0);
    std::vector<int> intervalCounts(nbStates, 0);
    arma::Mat<int> outCounts(nbStates, nbStates, arma::fill::zeros);

    int currentState = states[0];
    for (int i = 1; i < states.size(); ++i)
    {
      timeInStates[states[i - 1] - 1] += times[i] - times[i - 1];
      outCounts(states[i - 1] - 1, states[i] - 1)++;
      if (states[i] != currentState)
      {
        intervalCounts[currentState - 1]++;
        currentState = states[i];
      }
    }
    intervalCounts[currentState - 1]++; // Count the last interval

    // gibbs sampling
    std::vector<double> shape;
    std::vector<double> scale;
    shape.reserve(nbStates);
    scale.reserve(nbStates);
    std::transform(intervalCounts.begin(), intervalCounts.end(), priorShape.begin(),
                   std::back_inserter(shape), std::plus<double>());
    std::transform(timeInStates.begin(), timeInStates.end(), priorRate.begin(),
                   std::back_inserter(scale), std::plus<double>());
    std::transform(scale.begin(), scale.end(), scale.begin(), [](double x)
                   { return 1 / x; });

    for (int j = 0; j < nbStates; ++j) {
      // sample rates out of each state
      Q(j, j) = -std::min(Rcpp::rgamma(1, shape[j], scale[j])[0], kappa);
    }

    for (int j = 0; j < nbStates; ++j)
    {
      // sample probabilities into the other states (dirichlet random number generation)
      std::vector<double> probs(nbStates, 0);
      for (int k = 0; k < nbStates; ++k)
      {
        if (j != k)
        {
          double alpha = outCounts(j, k) + priorCon[k];
          probs[k] = Rcpp::rgamma(1, alpha, 1)[0];
        }
      }
      for (int k = 0; k < nbStates; ++k)
      {
        double sum = std::accumulate(probs.begin(), probs.end(), 0.0, [j, k](double acc, double val) { return k != j ? acc + val : acc; });
        if (j != k)
        {
          Q(j, k) = -Q(j, j) * probs[k] / sum;
        }
      }
    }
  }
  Rcpp::NumericMatrix getQ(const double time, const double lat, const double lng) const override
  {
    return Rcpp::wrap(Q);
  }
};

// Spatially aware NA Model
// In this transition rate model, individuals have their own rate generator matrices
// derived from Gibbs sampling on the number of intervals and the time spent in each state
// from their current state sequence.
// Transition probabilities into a ranged state are concentrated towards the nearest range center
// Transition probabilities are otherwise concentrated towards the nearest range center or unranged state.  
class SpatialNAModel : public Model
{
private:
  int nbStates;
  //NB times and states of switches only! not full data
  std::vector<double> times;
  std::vector<int> states;
  arma::vec mu;
  Rcpp::NumericVector priorShape;
  Rcpp::NumericVector priorRate;
  Rcpp::NumericVector priorCon;
  double kappa;
  arma::mat Q;

  double euclideanDistance(const double x1, const double y1, const double x2, const double y2) const {
    return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
  }

public:
  SpatialNAModel(const int nbStates, const Rcpp::DataFrame data, const arma::vec mu, const Rcpp::NumericVector priorShape, const Rcpp::NumericVector priorRate, const Rcpp::NumericVector priorCon, const double kappa)
      : nbStates(nbStates),
        times(Rcpp::as<std::vector<double>>(data["time"])),
        states(Rcpp::as<std::vector<int>>(data["state"])),
        mu(mu),
        priorShape(priorShape),
        priorRate(priorRate),
        priorCon(priorCon),
        kappa(kappa),
        Q(nbStates, nbStates, arma::fill::zeros)
  {
    std::vector<double> timeInStates(nbStates, 0.0);
    std::vector<int> intervalCounts(nbStates, 0);
    arma::Mat<int> outCounts(nbStates, nbStates, arma::fill::zeros);

    int currentState = states[0];
    for (int i = 1; i < states.size(); ++i)
    {
      timeInStates[states[i - 1] - 1] += times[i] - times[i - 1];
      outCounts(states[i - 1] - 1, states[i] - 1)++;
      if (states[i] != currentState)
      {
        intervalCounts[currentState - 1]++;
        currentState = states[i];
      }
    }
    intervalCounts[currentState - 1]++; // Count the last interval

    // gibbs sampling (rates only)
    std::vector<double> shape;
    std::vector<double> scale;
    shape.reserve(nbStates);
    scale.reserve(nbStates);
    std::transform(intervalCounts.begin(), intervalCounts.end(), priorShape.begin(),
                   std::back_inserter(shape), std::plus<double>());
    std::transform(timeInStates.begin(), timeInStates.end(), priorRate.begin(),
                   std::back_inserter(scale), std::plus<double>());
    std::transform(scale.begin(), scale.end(), scale.begin(), [](double x)
                   { return 1 / x; });

    for (int j = 0; j < nbStates; ++j) {
      // sample rates out of each state
      Q(j, j) = -std::min(Rcpp::rgamma(1, shape[j], scale[j])[0], kappa);
    }
  }
  Rcpp::NumericMatrix getQ(const double time, const double lat, const double lng) const override
  {
    // TODO: consider number of transitions (outCounts) in vicinity of current point?
    arma::mat Q_ = Q;
    arma::vec rangeDists(nbStates);

    for (int i = 0; i < nbStates; ++i) {
      double x = mu[(i * 2)];
      double y = mu[(i * 2) + 1];
      if(std::isfinite(x) && std::isfinite(y)) {
        rangeDists[i] = euclideanDistance(lng, lat, x, y);
      } else {
        rangeDists[i] = arma::datum::nan;
      }
    }
    double maxDist = arma::max(rangeDists);
    rangeDists.transform([maxDist](double dist) { return (std::isnan(dist) ? maxDist - maxDist : maxDist - dist); });

    for (int j = 0; j < nbStates; ++j)
    {
      // sample probabilities into the other states (dirichlet random number generation)
      std::vector<double> probs(nbStates, 0);
      for (int k = 0; k < nbStates; ++k)
      {
        if (j != k)
        {
          double alpha = rangeDists[k] + priorCon[k];
          probs[k] = Rcpp::rgamma(1, alpha, 1)[0];
        }
      }
      for (int k = 0; k < nbStates; ++k)
      {
        double sum = std::accumulate(probs.begin(), probs.end(), 0.0, [j, k](double acc, double val) { return k != j ? acc + val : acc; });
        if (j != k)
        {
          Q_(j, k) = -Q(j, j) * probs[k] / sum;
        }
      }
    }
    return Rcpp::wrap(Q_);
  }
};

// SST Model
class SSTModel : public Model
{
private:
  double sst_alpha, alpha;

public:
  // TODO: make these accept Named NumericVector
  SSTModel(double sst_alpha, double alpha) : sst_alpha(sst_alpha), alpha(alpha) {}
  Rcpp::NumericMatrix getQ(const double time, const double lat, const double lng) const override
  {
    return sst_alpha * std::exp(-alpha * time);
  }
};

// Time Model
class TimeModel : public Model
{
private:
  double t_alpha, alpha;

public:
  TimeModel(double t_alpha, double alpha) : t_alpha(t_alpha), alpha(alpha) {}
  Rcpp::NumericMatrix getQ(const double time, const double lat, const double lng) const override
  {
    return t_alpha * (1 - std::exp(-alpha * time));
  }
};

// Factory function to create the appropriate model
std::unique_ptr<Model> createModel(const int &nbStates, const std::string &model_type, const Rcpp::List &params, const double &kappa)
{
  if (model_type == "Null")
  {
    return std::make_unique<NullModel>(nbStates, params[0]);
  }
  else if (model_type == "NA")
  {
    Rcpp::DataFrame data(params["data"]);
    Rcpp::NumericVector priorShape(params["priorShape"]);
    Rcpp::NumericVector priorRate(params["priorRate"]);
    Rcpp::NumericVector priorCon(params["priorCon"]);
    return std::make_unique<NAModel>(nbStates, data, priorShape, priorRate, priorCon, kappa);
  }
  else if (model_type == "SpatialNA")
  {
    Rcpp::DataFrame data(params["data"]);
    arma::vec mu(Rcpp::as<arma::vec>(params["mu"]));
    Rcpp::NumericVector priorShape(params["priorShape"]);
    Rcpp::NumericVector priorRate(params["priorRate"]);
    Rcpp::NumericVector priorCon(params["priorCon"]);
    return std::make_unique<SpatialNAModel>(nbStates, data, mu, priorShape, priorRate, priorCon, kappa);
  }
  else
  {
    Rcpp::stop("Unknown model type");
  }
}