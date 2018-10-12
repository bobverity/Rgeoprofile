
#include "MCMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor for MCMC class
MCMC::MCMC() {
  
  double GTI_pow = 1.0; // TODO - set GTI power
  
  // thermodynamic parameters. The object beta_raised_vec stores values of beta
  // (the thermodynamic power), raised to the power GTI_pow
  beta_raised_vec = vector<double>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    beta_raised_vec[rung] = (rungs==1) ? 1 : pow((rung+1)/double(rungs), GTI_pow);
  }
  rung_order = seq_int(0,rungs-1);
  cold_rung = rung_order[rungs-1];
  
  // vector of particles
  particle_vec = vector<Particle>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    particle_vec[rung] = Particle(beta_raised_vec[rung]);
  }
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  
  // qmatrices
  log_qmatrix_running = vector<vector<double>>(n, vector<double>(K));
  if (K == 1) {
    qmatrix_final = vector<vector<double>>(n, vector<double>(K, samples));
  } else {
    qmatrix_final = vector<vector<double>>(n, vector<double>(K));
  }
  
  // objects for storing results
  loglike_burnin = vector<vector<double>>(rungs, vector<double>(burnin));
  loglike_sampling = vector<vector<double>>(rungs, vector<double>(samples));
  source_lon = vector<vector<double>>(samples, vector<double>(K));
  source_lat = vector<vector<double>>(samples, vector<double>(K));
  sigma = vector<vector<double>>(samples, vector<double>(K));
  expected_popsize = vector<double>(samples);
  
  // objects for storing acceptance rates
  source_accept = vector<int>(K);
  sigma_accept = vector<int>(K);
  coupling_accept = vector<int>(rungs-1);
  
  // store convergence
  rung_converged = vector<bool>(rungs, false);
  
}

//------------------------------------------------
// run burn-in phase of MCMC
void MCMC::burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // print header
  if (!silent) {
    print("Running MCMC for K =", K);
    print("Burn-in phase");
  }
  
  // read in R functions
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // define points at which convergence checked
  vector<int> convergence_checkpoint(1,converge_test);
  while(convergence_checkpoint.back()<burnin) {
    convergence_checkpoint.push_back(convergence_checkpoint.back()+converge_test);
  }
  int checkpoint_i = 0;
  
  // reset particles
  for (int r=0; r<rungs; r++) {
    particle_vec[r].reset();
    particle_vec[r].beta_raised = beta_raised_vec[r];
  }
  rung_order = seq_int(0,rungs-1);
  
  // loop through burnin iterations
  vector<bool> convergence_reached(rungs, false);
  bool all_convergence_reached = false;
  for (int rep=0; rep<burnin; rep++) {
    
    // update particles
    for (int r=0; r<rungs; r++) {
      
      // skip over converged rungs
      if (convergence_reached[r]) {
        continue;
      }
      int rung = rung_order[r];
      
      // update sources
      particle_vec[rung].update_sources(true, rep+1);
      
      // update sigma
      particle_vec[rung].update_sigma(true, rep+1);
      
      // update expected population size
      particle_vec[rung].update_expected_popsize();
      
    } // end loop over rungs
    
    // focus on cold rung
    cold_rung = rung_order[rungs-1];
    
    // methods that only apply when K>1
    if (K > 1) {
      
      // update qmatrix of cold rung
      particle_vec[cold_rung].update_qmatrix();
      
      // fix labels
      particle_vec[cold_rung].solve_label_switching(log_qmatrix_running);
      label_order = particle_vec[cold_rung].label_order;
      
      // add particle log_qmatrix to log_qmatrix_running
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          log_qmatrix_running[i][k] = log_sum(log_qmatrix_running[i][k], particle_vec[cold_rung].log_qmatrix[i][label_order[k]]);
        }
      }
      
    }
    
    // store loglikelihood
    for (int r=0; r<rungs; r++) {
      if (convergence_reached[r]) {
        continue;
      }
      int rung = rung_order[r];
      loglike_burnin[r][rep] = particle_vec[rung].loglike;
    }
    
    // update progress bars
    if (!silent) {
      if ((rep+1) == burnin) {
        update_progress(args_progress, "pb_burnin", rep+1, burnin);
      } else {
        int remainder = rep % int(ceil(double(burnin)/100));
        if (remainder == 0 && !pb_markdown) {
          update_progress(args_progress, "pb_burnin", rep+1, burnin);
        }
      }
    }
    
    // check for convergence
    if (auto_converge && (rep+1) == convergence_checkpoint[checkpoint_i]) {
      
      // check for convergence of all unconverged chains
      for (int r=0; r<rungs; r++) {
        if (!convergence_reached[r]) {
          convergence_reached[r] = rcpp_to_bool(test_convergence(loglike_burnin[r], rep+1));
          if (convergence_reached[r]) {
            rung_converged[r] = true;
            loglike_burnin[r].resize(rep+1);
          }
        }
      }
      // break if convergence reached
      all_convergence_reached = true;
      for (int r=0; r<rungs; r++) {
        if (!convergence_reached[r]) {
          all_convergence_reached = false;
          break;
        }
      }
      // end if all reached convergence
      if (all_convergence_reached) {
        if (!silent) {
          update_progress(args_progress, "pb_burnin", burnin, burnin);
          print("   converged within", rep+1, "iterations");
        }
        break;
      }
      checkpoint_i++;
      
    }  // end if auto_converge
    
  }  // end burn-in iterations
  
  // warning if still not converged
  if (!all_convergence_reached && !silent) {
    print("   Warning: convergence still not reached within", burnin, "iterations");
  }
  
}

//------------------------------------------------
// run sampling phase of MCMC
void MCMC::sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // print header
  if (!silent) {
    print("Sampling phase");
  }
  
  // read in R functions
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // loop through sampling iterations
  for (int rep=0; rep<samples; rep++) {
    
    // update particles
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // update sources
      particle_vec[rung].update_sources(false, 0);
      
      // update sigma
      particle_vec[rung].update_sigma(false, 0);
      
      // update expected population size
      particle_vec[rung].update_expected_popsize();
      
    } // end loop over rungs
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // focus on cold rung
    cold_rung = rung_order[rungs-1];
    
    // methods that only apply when K>1
    if (K > 1) {
      
      // update qmatrix of cold rung
      particle_vec[cold_rung].update_qmatrix();
      
      // fix labels
      particle_vec[cold_rung].solve_label_switching(log_qmatrix_running);
      label_order = particle_vec[cold_rung].label_order;
      
      // add particle log_qmatrix to log_qmatrix_running
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          log_qmatrix_running[i][k] = log_sum(log_qmatrix_running[i][k], particle_vec[cold_rung].log_qmatrix[i][label_order[k]]);
        }
      }
      
      // add particle qmatrix to qmatrix_final
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          qmatrix_final[i][k] += particle_vec[cold_rung].qmatrix[i][label_order[k]];
        }
      }
      
    }
    
    // store loglikelihood
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      loglike_sampling[r][rep] = particle_vec[rung].loglike;
    }
    
    // store source locations
    for (int k=0; k<K; k++) {
      source_lon[rep][k] = particle_vec[cold_rung].source_lon[label_order[k]];
      source_lat[rep][k] = particle_vec[cold_rung].source_lat[label_order[k]];
    }
    
    // store sigma
    for (int k=0; k<K; k++) {
      sigma[rep][k] = particle_vec[cold_rung].sigma[label_order[k]];
    }
    
    // store expected population size
    expected_popsize[rep] = particle_vec[cold_rung].expected_popsize;
    
    // update progress bars
    if (!silent) {
      if ((rep+1) == samples) {
        update_progress(args_progress, "pb_samples", rep+1, samples);
      } else {
        int remainder = rep % int(ceil(double(samples)/100));
        if (remainder == 0 && !pb_markdown) {
          update_progress(args_progress, "pb_samples", rep+1, samples);
        }
      }
    }
    
  } // end sampling iterations
  
  // store acceptance rates
  source_accept = particle_vec[cold_rung].source_accept;
  sigma_accept = particle_vec[cold_rung].sigma_accept;
  
}

//------------------------------------------------
// Metropolis-coupling to propose swaps between temperature rungs
void MCMC::metropolis_coupling() {
  
  // loop over rungs, starting with the hottest chain and moving to the cold
  // chain. Each time propose a swap with the next rung up
  for (int i=0; i<(rungs-1); i++) {
    
    // define rungs of interest
    int rung1 = rung_order[i];
    int rung2 = rung_order[i+1];
    
    // get log-likelihoods and beta values of two chains in the comparison
    double loglike1 = particle_vec[rung1].loglike;
    double loglike2 = particle_vec[rung2].loglike;
    
    double beta_raised1 = particle_vec[rung1].beta_raised;
    double beta_raised2 = particle_vec[rung2].beta_raised;
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (loglike2*beta_raised1 + loglike1*beta_raised2) - (loglike1*beta_raised1 + loglike2*beta_raised2);

    // accept or reject move
    double rand1 = runif1();
    if (log(rand1)<acceptance) {
      
      // swap beta values
      particle_vec[rung1].beta_raised = beta_raised2;
      particle_vec[rung2].beta_raised = beta_raised1;
      
      // TODO - swap proposal SDs
      
      // swap rung order
      rung_order[i] = rung2;
      rung_order[i+1] = rung1;
      
      // update acceptance rates
      coupling_accept[i]++;
    }
  }
  
}
