
#include <chrono>

#include "main.h"
#include "Data.h"
#include "Parameters.h"
#include "Spatial_prior.h"
#include "Lookup.h"
#include "MCMC.h"
#include "misc.h"
//#include "probability.h"
//#include "hungarian.h"

using namespace std;

//------------------------------------------------
// run main MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_data = args["args_data"];
  Rcpp::List args_model = args["args_model"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  
  // read in data into separate class
  Data data(args_data);
  
  // read in parameters into separate class
  Parameters parameters(args_model);
  
  // define spatial prior
  Spatial_prior spatial_prior(args_model);
  /*
  vector<vector<double>> z(parameters.n_lat, vector<double>(parameters.n_lon));
  int i2 = 0;
  for (int i=0; i<parameters.n_lat; ++i) {
    for (int j=0; j<parameters.n_lon; ++j) {
      //z[i][j] = spatial_prior.spatial_prior_mask[i2++];
    }
  }
  
  Rcpp::List ret2 = Rcpp::List::create(Rcpp::Named("z") = z);
  return(ret2);
  */
  // define lookup table from data and parameters. The lookup table 
  // precalculates the distance in km of every data point from every possible 
  // source location in a lon/lat grid. The grid is defined from the lon/lat min
  // and max values ()in the parameter set), and the lon/lat precision (in the
  // MCMC arguments).
  //Lookup lookup(data, parameters);
  Lookup lookup;
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create MCMC object
  MCMC mcmc;
  
  // run MCMC
  mcmc.burnin_mcmc(args_functions, args_progress);
  mcmc.sampling_mcmc(args_functions, args_progress);
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!parameters.silent) {
    print("   completed in", time_span.count(), "seconds\n");
  }
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( mcmc.loglike_burnin ));
  ret.push_back(Rcpp::wrap( mcmc.loglike_sampling ));
  ret.push_back(Rcpp::wrap( mcmc.source_lon ));
  ret.push_back(Rcpp::wrap( mcmc.source_lat ));
  ret.push_back(Rcpp::wrap( mcmc.sigma ));
  ret.push_back(Rcpp::wrap( mcmc.expected_popsize ));
  ret.push_back(Rcpp::wrap( mcmc.qmatrix_final ));
  ret.push_back(Rcpp::wrap( mcmc.source_accept ));
  ret.push_back(Rcpp::wrap( mcmc.sigma_accept ));
  ret.push_back(Rcpp::wrap( mcmc.rung_converged ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("loglike_burnin");
  ret_names.push_back("loglike_sampling");
  ret_names.push_back("source_lon");
  ret_names.push_back("source_lat");
  ret_names.push_back("sigma");
  ret_names.push_back("expected_popsize");
  ret_names.push_back("qmatrix");
  ret_names.push_back("source_accept");
  ret_names.push_back("sigma_accept");
  ret_names.push_back("rung_converged");
  
  ret.names() = ret_names;
  return ret;
}
