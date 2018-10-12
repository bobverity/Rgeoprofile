
#include "Parameters.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// declare static member variables

// MCMC parameters
int Parameters::burnin;
int Parameters::samples;
int Parameters::rungs;
bool Parameters::auto_converge;
int Parameters::converge_test;
bool Parameters::pb_markdown;
bool Parameters::silent;
bool Parameters::coupling_on;

// model parameters
int Parameters::K;
double Parameters::sentinel_radius;
double Parameters::min_lon;
double Parameters::max_lon;
double Parameters::res_lon;
int Parameters::n_lon;
double Parameters::min_lat;
double Parameters::max_lat;
double Parameters::res_lat;
int Parameters::n_lat;
vector<double> Parameters::source_init;
int Parameters::sigma_model;
double Parameters::sigma_prior_meanlog;
double Parameters::sigma_prior_sdlog;
double Parameters::expected_popsize_prior_sd;
double Parameters::expected_popsize_prior_shape;
double Parameters::expected_popsize_prior_rate;

int Parameters::index_lon_min;
int Parameters::index_lat_min;

//------------------------------------------------
// constructor for Parameters class
Parameters::Parameters(const Rcpp::List &args) {
  
  // MCMC parameters
  burnin = rcpp_to_int(args["burnin"]);
  samples = rcpp_to_int(args["samples"]);
  rungs = 1;  // TODO - variable number of rungs
  auto_converge = rcpp_to_bool(args["auto_converge"]);
  converge_test = rcpp_to_int(args["converge_test"]);
  pb_markdown = rcpp_to_bool(args["pb_markdown"]);
  silent = rcpp_to_bool(args["silent"]);
  coupling_on = false;  // TODO - metropolis-coupling over temperature rungs
  
  // model parameters
  K = rcpp_to_int(args["K"]);
  sentinel_radius = rcpp_to_double(args["sentinel_radius"]);
  min_lon = rcpp_to_double(args["min_lon"]);
  max_lon = rcpp_to_double(args["max_lon"]);
  res_lon = rcpp_to_double(args["res_lon"]);
  n_lon = rcpp_to_int(args["n_lon"]);
  min_lat = rcpp_to_double(args["min_lat"]);
  max_lat = rcpp_to_double(args["max_lat"]);
  res_lat = rcpp_to_double(args["res_lat"]);
  n_lat = rcpp_to_int(args["n_lat"]);
  source_init = rcpp_to_vector_double(args["source_init"]);
  sigma_model = rcpp_to_int(args["sigma_model_numeric"]);
  
  // get sigma prior mean and sd in log space from raw inputs
  double sigma_prior_mean = rcpp_to_double(args["sigma_prior_mean"]);
  double sigma_prior_sd = rcpp_to_double(args["sigma_prior_sd"]);
  double sigma_prior_varlog = log(pow(sigma_prior_sd,2)/pow(sigma_prior_mean,2) + 1);
  sigma_prior_sdlog = sqrt(sigma_prior_varlog);
  sigma_prior_meanlog = log(sigma_prior_mean) - sigma_prior_varlog/2.0;
  
  // get expected_popsize shape and rate parameters from raw inputs
  double expected_popsize_prior_mean = rcpp_to_double(args["expected_popsize_prior_mean"]);
  expected_popsize_prior_sd = rcpp_to_double(args["expected_popsize_prior_sd"]);
  if (expected_popsize_prior_sd == -1) {  // improper prior
    expected_popsize_prior_shape = 0;
    expected_popsize_prior_rate = 0;
  } else {  // proper prior
    expected_popsize_prior_shape = pow(expected_popsize_prior_mean,2)/pow(expected_popsize_prior_sd,2);
    expected_popsize_prior_rate = expected_popsize_prior_mean/pow(expected_popsize_prior_sd,2);
  }
  
  // get index limits. These are used when converting continuous lon/lat values
  // to an integer index starting at zero
  index_lon_min = min_lon/res_lon;
  index_lat_min = min_lat/res_lat;
  
}
