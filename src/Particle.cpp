
#include "Particle.h"
#include "misc.h"
#include "probability.h"
#include "hungarian.h"

using namespace std;

//------------------------------------------------
// constructor for Particle class
Particle::Particle(double beta_raised) {
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  this->beta_raised = beta_raised;
  
  // source locations
  source_lon = vector<double>(K);
  source_lat = vector<double>(K);
  
  // standard deviation of sources (km)
  sigma = vector<double>(K, 1.0);
  
  // scaling factor on hazard surface, equivalent to the expected total
  // population size (both observed and unobserved) in unit time
  expected_popsize = 500.0;
  log_expected_popsize = log(expected_popsize);
  
  // qmatrices
  log_qmatrix = vector<vector<double>>(n, vector<double>(K));
  qmatrix = vector<vector<double>>(n, vector<double>(K));
  
  // proposal standard deviations
  source_propSD = vector<double>(K, 1.0);
  sigma_propSD = vector<double>(K, 1.0);
  
  // Robbins-Monro stepsize
  source_rm_stepsize = 1.0;
  sigma_rm_stepsize = 1.0;
  
  // misc constants
  // area around a sentinel site
  log_sentinel_area = LOG_PI + 2*log(sentinel_radius);
  // area of source prior
  log_search_area = log(max_lon-min_lon) + log(max_lat-min_lat);
  // sum of counts over all sentinel sites
  counts_total = sum(sentinel_counts);
  // log of K
  log_K = log(K);
  
  // likelihood
  dist_source_data = vector<vector<double>>(n, vector<double>(K));
  dist_source_data_prop = vector<double>(n);
  log_hazard_height = vector<vector<double>>(n, vector<double>(K));
  log_hazard_height_prop = vector<double>(n);
  log_hazard_height_prop2 = vector<vector<double>>(n, vector<double>(K));
  loglike = 0;
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  label_order_new = vector<int>(K);
  
  // objects for solving label switching problem
  cost_mat = vector<vector<double>>(K, vector<double>(K));
  best_perm = vector<int>(K);
  best_perm_order = vector<int>(K);
  edges_left = vector<int>(K);
  edges_right = vector<int>(K);
  blocked_left = vector<int>(K);
  blocked_right = vector<int>(K);
  
  // store acceptance rates
  source_accept = vector<int>(K);
  sigma_accept = vector<int>(K);
  
}

//------------------------------------------------
// reset particle
void Particle::reset() {
  
  // initialise source locations
  for (int k=0; k<K; ++k) {
    source_lon[k] = source_init[0];
    source_lat[k] = source_init[1];
  }
  
  // draw sigma from prior
  for (int k=0; k<K; ++k) {
    sigma[k] = 1.0; // TODO - draw from prior
  }
  
  // draw expected popsize from prior
  expected_popsize = 100.0;  // TODO - draw from prior
  log_expected_popsize = log(expected_popsize);
  
  // initialise proposal standard deviations
  source_propSD = vector<double>(K, 0.01);
  sigma_propSD = vector<double>(K, 1.0);
  
  // calculate initial likelihood. Calling calculate_loglike_source() on each 
  // source updates the dist_source_data_prop and log_hazard_height_prop 
  // objects, which can then be stored back into the final matrices. This is 
  // equivalent to running a Metropolis-Hastings step in which the move is 
  // guaranteed to be accepted
  for (int k=0; k<K; ++k) {
    loglike = calculate_loglike_source(source_lon[k], source_lat[k], k);
    for (int i=0; i<n; ++i) {
      dist_source_data[i][k] = dist_source_data_prop[i];
      log_hazard_height[i][k] = log_hazard_height_prop[i];
    }
  }
  
  // reset acceptance rates
  source_accept = vector<int>(K);
  sigma_accept = vector<int>(K);
  
}

//------------------------------------------------
// calculate log-likelihood given new proposed source
double Particle::calculate_loglike_source(double source_lon_prop, double source_lat_prop, int k) {
  
  // get prior probability
  double prior_prob = get_value(source_lon_prop, source_lat_prop);
  if (prior_prob == 0) {
    return -OVERFLO;
  }
  
  // initialise new likelihood
  double loglike_prop = log(prior_prob);
  
  // loop through sentinel sites
  for (int i=0; i<n; ++i) {
    
    // get distance from proposed source to data point i
    double dist = get_data_dist(source_lon_prop, source_lat_prop, i);
    dist_source_data_prop[i] = dist;
    
    // calculate bivariate normal height of data point i from proposed source. 
    // This is equivalent to the univariate normal density of the distance 
    // between data and source in lon, multiplied by the univariate normal 
    // density of the distance between data and source in lat. Due to the 
    // properties of normal distributions this is equivalent to the univariate 
    // normal density of the euclidian distance between source and data, 
    // multiplied by the univariate normal density of 0 from 0 (the latter is 
    // needed to make it a bivariate density). Finally, in log space densities 
    // are summed not multiplied.
    log_hazard_height_prop[i] = dnorm1(dist, 0, sigma[k], true) + dnorm1(0, 0, sigma[k], true);
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop[i];
    for (int j=0; j<K; ++j) {
      if (j == k) {
        continue;
      }
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // divide by K
    log_hazard_sum -= log_K;
    
    // define the rate lambda of the Poisson process at this sentinel site,
    // while remaining in log space
    double log_lambda = log_sentinel_area + log_expected_popsize + log_hazard_sum;
    
    // calculate the Poisson log-probability of the counts at this sentinel site
    loglike_prop += sentinel_counts[i]*log_lambda - exp(log_lambda) - lgamma(sentinel_counts[i]+1);
  }
  
  return loglike_prop;
}

//------------------------------------------------
// update source locations
void Particle::update_sources(bool robbins_monro_on, int iteration) {
  
  // loop through all sources
  for (int k=0; k<K; ++k) {
    
    // propose new source location
    double source_lon_prop = rnorm1(source_lon[k], source_propSD[k]);
    double source_lat_prop = rnorm1(source_lat[k], source_propSD[k]);
    
    // check proposed source within defined range
    if (source_lon_prop <= min_lon || source_lon_prop >= max_lon || source_lat_prop <= min_lat || source_lat_prop >= max_lat) {
      
      // auto-reject proposed move
      if (robbins_monro_on) {
        source_propSD[k] = exp(log(source_propSD[k]) - source_rm_stepsize*0.234/sqrt(iteration));
      }
      continue;
    }
    
    // calculate new loglikelihood
    double loglike_prop = calculate_loglike_source(source_lon_prop, source_lat_prop, k);
    
    // Metropolis-Hastings ratio
    double MH_ratio = loglike_prop - loglike;
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update source
      source_lon[k] = source_lon_prop;
      source_lat[k] = source_lat_prop;
      
      // update stored distances and hazard values
      for (int i=0; i<n; ++i) {
        dist_source_data[i][k] = dist_source_data_prop[i];
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        source_propSD[k] = exp(log(source_propSD[k]) + sigma_rm_stepsize*(1-0.23)/sqrt(iteration));
      } else {
        source_accept[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        source_propSD[k] = exp(log(source_propSD[k]) - sigma_rm_stepsize*0.23/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
  }  // end loop through sources
  
}

//------------------------------------------------
// update sigma
void Particle::update_sigma(bool robbins_monro_on, int iteration) {
  
  // update single sigma or separately for each source
  if (sigma_model == 1) {
    update_sigma_single(robbins_monro_on, iteration);
  } else if (sigma_model == 2) {
    update_sigma_independent(robbins_monro_on, iteration);
  }
  
}

//------------------------------------------------
// update single sigma for all sources
void Particle::update_sigma_single(bool robbins_monro_on, int iteration) {
  
  // propose new value
  double sigma_prop = rnorm1(sigma[0], sigma_propSD[0]);
  sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
  
  // initialise new likelihood
  double loglike_prop = 0;
  
  // loop through sentinel sites
  for (int i=0; i<n; ++i) {
    
    // loop through sources
    for (int k=0; k<K; ++k) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop2[i][k] = dnorm1(dist, 0, sigma_prop, true) + dnorm1(0, 0, sigma_prop, true);
    }
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = log_hazard_height_prop2[i][0];
    for (int j=1; j<K; ++j) {
      if (log_hazard_sum < log_hazard_height_prop2[i][j]) {
        log_hazard_sum = log_hazard_height_prop2[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height_prop2[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height_prop2[i][j] - log_hazard_sum));
      }
    }
    
    // divide by K
    log_hazard_sum -= log_K;
    
    // define the rate lambda of the Poisson process at this sentinel site,
    // while remaining in log space
    double log_lambda = log_sentinel_area + log_expected_popsize + log_hazard_sum;
    
    // calculate the Poisson log-probability of the counts at this sentinel site
    loglike_prop += sentinel_counts[i]*log_lambda - exp(log_lambda) - lgamma(sentinel_counts[i]+1);
    
  }
  
  // calculate priors
  double logprior = dlnorm1(sigma[0], sigma_prior_meanlog, sigma_prior_sdlog);
  double logprior_prop = dlnorm1(sigma_prop, sigma_prior_meanlog, sigma_prior_sdlog);
  
  // Metropolis-Hastings ratio
  double MH_ratio = (loglike_prop + logprior_prop) - (loglike + logprior);
  
  // Metropolis-Hastings step
  if (log(runif_0_1())<MH_ratio) {
    
    // update sigma for all sources
    for (int k=0; k<K; ++k) {
      sigma[k] = sigma_prop;
    }
    
    // update stored hazard values
    log_hazard_height = log_hazard_height_prop2;
    
    // update likelihood
    loglike = loglike_prop;
    
    // Robbins-Monro positive update (on the log scale)
    if (robbins_monro_on) {
      sigma_propSD[0] = exp(log(sigma_propSD[0]) + sigma_rm_stepsize*(1-0.23)/sqrt(iteration));
    } else {
      sigma_accept[0]++;
    }
    
  } else {
    
    // Robbins-Monro negative update (on the log scale)
    if (robbins_monro_on) {
      sigma_propSD[0] = exp(log(sigma_propSD[0]) - sigma_rm_stepsize*0.23/sqrt(iteration));
    }
    
  }  // end Metropolis-Hastings step
  
}

//------------------------------------------------
// update independent sigma for each source
void Particle::update_sigma_independent(bool robbins_monro_on, int iteration) {
  
  // return if prior is exact
  if (sigma_prior_sdlog == 0) {
    return;
  }
  
  // loop through sources
  for (int k=0; k<K; ++k) {
    
    // propose new value
    double sigma_prop = rnorm1(sigma[k], sigma_propSD[k]);
    sigma_prop = (sigma_prop < UNDERFLO) ? UNDERFLO : sigma_prop;
    
    // initialise new likelihood
    double loglike_prop = 0;
    
    // loop through sentinel sites
    for (int i=0; i<n; ++i) {
      
      // recalculate hazard given new sigma
      double dist = dist_source_data[i][k];
      log_hazard_height_prop[i] = dnorm1(dist, 0, sigma_prop, true) + dnorm1(0, 0, sigma_prop, true);
      
      // sum hazard over sources while remaining in log space
      double log_hazard_sum = log_hazard_height_prop[i];
      for (int j=0; j<K; ++j) {
        if (j == k) {
          continue;
        }
        if (log_hazard_sum < log_hazard_height[i][j]) {
          log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
        } else {
          log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
        }
      }
      
      // divide by K
      log_hazard_sum -= log_K;
      
      // define the rate lambda of the Poisson process at this sentinel site,
      // while remaining in log space
      double log_lambda = log_sentinel_area + log_expected_popsize + log_hazard_sum;
      
      // calculate the Poisson log-probability of the counts at this sentinel site
      loglike_prop += sentinel_counts[i]*log_lambda - exp(log_lambda) - lgamma(sentinel_counts[i]+1);
    }
    
    // calculate priors
    double logprior = dlnorm1(sigma[k], sigma_prior_meanlog, sigma_prior_sdlog);
    double logprior_prop = dlnorm1(sigma_prop, sigma_prior_meanlog, sigma_prior_sdlog);
    
    // Metropolis-Hastings ratio
    double MH_ratio = (loglike_prop + logprior_prop) - (loglike + logprior);
    
    // Metropolis-Hastings step
    if (log(runif_0_1()) < MH_ratio) {
      
      // update sigma for this source
      sigma[k] = sigma_prop;
      
      // update stored hazard values
      for (int i=0; i<n; ++i) {
        log_hazard_height[i][k] = log_hazard_height_prop[i];
      }
      
      // update likelihood
      loglike = loglike_prop;
      
      // Robbins-Monro positive update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) + sigma_rm_stepsize*(1-0.234)/sqrt(iteration));
      } else {
        sigma_accept[k]++;
      }
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (robbins_monro_on) {
        sigma_propSD[k] = exp(log(sigma_propSD[k]) - sigma_rm_stepsize*0.234/sqrt(iteration));
      }
      
    }  // end Metropolis-Hastings step
    
  }  // end loop over sources
  
}

//------------------------------------------------
// update separate sigma for each source
void Particle::update_expected_popsize() {
  
  // return if prior is exact
  if (expected_popsize_prior_sd == 0) {
    return;
  }
  
  // sum of Poisson rate over sentinel sites
  double lambda_total = 0;
  for (int i=0; i<n; ++i) {
    
    // take mean of hazard over sources
    for (int k=0; k<K; ++k) {
      lambda_total += exp(log_sentinel_area + log_hazard_height[i][k] - log_K);
    }
  }
  
  // draw new expected population size
  double posterior_shape = expected_popsize_prior_shape + counts_total;
  double posterior_rate = expected_popsize_prior_rate + lambda_total;
  expected_popsize = rgamma1(posterior_shape, posterior_rate);
  
}

//------------------------------------------------
// update qmatrix
void Particle::update_qmatrix() {
  
  // loop through sentinel sites
  for (int i=0; i<n; ++i) {
    
    // skip if no observations at this site
    if (sentinel_counts[i] == 0) {
      continue;
    }
    
    // sum hazard over sources while remaining in log space
    double log_hazard_sum = -OVERFLO;
    for (int j=0; j<K; ++j) {
      if (log_hazard_sum < log_hazard_height[i][j]) {
        log_hazard_sum = log_hazard_height[i][j] + log(1 + exp(log_hazard_sum - log_hazard_height[i][j]));
      } else {
        log_hazard_sum = log_hazard_sum + log(1 + exp(log_hazard_height[i][j] - log_hazard_sum));
      }
    }
    
    // qmatrix equal to normalised hazard
    for (int j=0; j<K; ++j) {
      log_qmatrix[i][j] = log_hazard_height[i][j] - log_hazard_sum;
      qmatrix[i][j] = exp(log_hazard_height[i][j] - log_hazard_sum);
    }
    
  }  // end loop through sentinel sites
  
}

//------------------------------------------------
// solve label switching problem
void Particle::solve_label_switching(const vector<vector<double>> &log_qmatrix_running) {
  
  // recalculate cost matrix
  for (int k1=0; k1<K; k1++) {
    fill(cost_mat[k1].begin(), cost_mat[k1].end(), 0);
    for (int k2=0; k2<K; k2++) {
      for (int i=0; i<n; i++) {
        if (sentinel_counts[i] > 0) {
          cost_mat[k1][k2] += qmatrix[i][label_order[k1]]*(log_qmatrix[i][label_order[k1]] - log_qmatrix_running[i][k2]);
        }
      }
    }
  }
  
  // find best permutation of current labels using Hungarian algorithm
  best_perm = hungarian(cost_mat, edges_left, edges_right, blocked_left, blocked_right);
  
  // define best_perm_order
  for (int k=0; k<K; k++) {
    best_perm_order[best_perm[k]] = k;
  }
  
  // replace old label order with new
  for (int k=0; k<K; k++) {
    label_order_new[k] = label_order[best_perm_order[k]];
  }
  label_order = label_order_new;
  
}

