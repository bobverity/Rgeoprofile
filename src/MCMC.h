
#pragma once

#include <Rcpp.h>
#include "Data.h"
#include "Parameters.h"
#include "Particle.h"

//------------------------------------------------
// class defining MCMC
class MCMC : public Data, public Parameters {
  
public:
  
  // thermodynamic parameters
  std::vector<double> beta_raised_vec;
  std::vector<int> rung_order;
  int cold_rung;
  
  // vector of particles
  std::vector<Particle> particle_vec;
  
  // ordering of labels
  std::vector<int> label_order;
  
  // qmatrices
  std::vector<std::vector<double>> log_qmatrix_running;
  std::vector<std::vector<double>> qmatrix_final;
  
  // objects for storing results
  std::vector<std::vector<double>> loglike_burnin;
  std::vector<std::vector<double>> loglike_sampling;
  std::vector<std::vector<double>> source_lon;
  std::vector<std::vector<double>> source_lat;
  std::vector<std::vector<double>> sigma;
  std::vector<double> expected_popsize;
  
  // objects for storing acceptance rates
  std::vector<int> source_accept;
  std::vector<int> sigma_accept;
  std::vector<int> coupling_accept;
  
  // store convergence
  std::vector<bool> rung_converged;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC();
  
  // other functions
  void burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void metropolis_coupling();
  
};
