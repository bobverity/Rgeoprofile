
#include "Parameters.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// declare static member variables

// MCMC parameters
double Parameters::precision_lon;
double Parameters::precision_lat;
int Parameters::burnin;
int Parameters::samples;
int Parameters::rungs;
bool Parameters::auto_converge;
int Parameters::converge_test;
bool Parameters::solve_label_switching_on;
bool Parameters::pb_markdown;
bool Parameters::silent;
bool Parameters::coupling_on;

// model parameters
int Parameters::K;
double Parameters::sentinel_radius;
double Parameters::source_min_lon;
double Parameters::source_max_lon;
double Parameters::source_min_lat;
double Parameters::source_max_lat;
int Parameters::sigma_model;
double Parameters::sigma_prior_meanlog;
double Parameters::sigma_prior_sdlog;
double Parameters::expected_popsize_prior_shape;
double Parameters::expected_popsize_prior_rate;

//------------------------------------------------
// constructor for Parameters class
Parameters::Parameters(const Rcpp::List &args) {
  
  // MCMC parameters
  precision_lon = rcpp_to_double(args["precision_lon"]);
  precision_lat = rcpp_to_double(args["precision_lat"]);
  burnin = rcpp_to_int(args["burnin"]);
  samples = rcpp_to_int(args["samples"]);
  rungs = 1;  // TODO - variable number of rungs
  auto_converge = rcpp_to_bool(args["auto_converge"]);
  converge_test = rcpp_to_int(args["converge_test"]);
  solve_label_switching_on = rcpp_to_bool(args["solve_label_switching_on"]);
  pb_markdown = rcpp_to_bool(args["pb_markdown"]);
  silent = rcpp_to_bool(args["silent"]);
  coupling_on = false;
  
  // model parameters
  K = rcpp_to_int(args["K"]);
  sentinel_radius = rcpp_to_double(args["sentinel_radius"]);
  source_min_lon = rcpp_to_double(args["source_min_lon"]);
  source_max_lon = rcpp_to_double(args["source_max_lon"]);
  source_min_lat = rcpp_to_double(args["source_min_lat"]);
  source_max_lat = rcpp_to_double(args["source_max_lat"]);
  sigma_model = rcpp_to_int(args["sigma_model_numeric"]);
  
  double sigma_prior_mean = rcpp_to_double(args["sigma_prior_mean"]);
  double sigma_prior_sd = rcpp_to_double(args["sigma_prior_sd"]);
  double sigma_prior_varlog = log(pow(sigma_prior_sd,2)/pow(sigma_prior_mean,2) + 1);
  sigma_prior_sdlog = sqrt(sigma_prior_varlog);
  sigma_prior_meanlog = log(sigma_prior_mean) - sigma_prior_varlog/2.0;
  
  expected_popsize_prior_shape = rcpp_to_double(args["expected_popsize_prior_shape"]);
  expected_popsize_prior_rate = rcpp_to_double(args["expected_popsize_prior_rate"]);
  
}
