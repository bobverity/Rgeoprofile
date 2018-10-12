
#pragma once

#include <Rcpp.h>

#include "Data.h"
#include "Parameters.h"

//------------------------------------------------
// class containing lookup tables
class Lookup : public Data, public Parameters {
  
public:
  
  // lookup table
  static std::vector<double> lookup_dist;
  
  // constructors
  Lookup();
  
  // other methods
  double get_data_dist(double source_lon, double source_lat, int data_index);
  
};
