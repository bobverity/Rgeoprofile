
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing all parameters
class Parameters {
  
public:
  
  // MCMC parameters
  static double precision_lon;
  static double precision_lat;
  static int burnin;
  static int samples;
  static int rungs;
  static bool auto_converge;
  static int converge_test;
  static bool solve_label_switching_on;
  static bool pb_markdown;
  static bool silent;
  static bool coupling_on;
  
  // model parameters
  static int K;
  static double sentinel_radius;
  static double min_lon;
  static double max_lon;
  static double min_lat;
  static double max_lat;
  static int sigma_model;
  static double sigma_prior_meanlog;
  static double sigma_prior_sdlog;
  static double expected_popsize_prior_sd;
  static double expected_popsize_prior_shape;
  static double expected_popsize_prior_rate;
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
};
