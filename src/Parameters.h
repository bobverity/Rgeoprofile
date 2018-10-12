
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing all parameters
class Parameters {
  
public:
  
  // MCMC parameters
  static int burnin;
  static int samples;
  static int rungs;
  static bool auto_converge;
  static int converge_test;
  static bool pb_markdown;
  static bool silent;
  static bool coupling_on;
  
  // model parameters
  static int K;
  static double sentinel_radius;
  static double min_lon;
  static double max_lon;
  static double res_lon;
  static int n_lon;
  static double min_lat;
  static double max_lat;
  static double res_lat;
  static int n_lat;
  static std::vector<double> source_init;
  static int sigma_model;
  static double sigma_prior_meanlog;
  static double sigma_prior_sdlog;
  static double expected_popsize_prior_sd;
  static double expected_popsize_prior_shape;
  static double expected_popsize_prior_rate;
  
  static int index_lon_min;
  static int index_lat_min;
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
};
