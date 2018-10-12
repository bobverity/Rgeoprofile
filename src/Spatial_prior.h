
#pragma once

#include <Rcpp.h>
#include "Parameters.h"

//------------------------------------------------
// class containing spatial prior
class Spatial_prior : public Parameters {
  
public:
  
  // spatial prior
  static std::vector<bool> spatial_prior_mask;
  
  // constructors
  Spatial_prior() {};
  Spatial_prior(const Rcpp::List &args);
  
  // other methods
  double get_value(double lon, double lat);
  
};
