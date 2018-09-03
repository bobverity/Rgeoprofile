
#pragma once

#include <Rcpp.h>

#include "Data.h"
#include "Parameters.h"

//------------------------------------------------
// class containing lookup tables
class Lookup {
  
public:
  
  // copy of some parameters
  static double precision_lon;
  static double precision_lat;
  static int n_data;
  
  // index limits
  static int index_limit_lon_min;
  static int index_limit_lon_max;
  static int index_limit_lat_min;
  static int index_limit_lat_max;
  static int n_lon;
  static int n_lat;
  
  // lookup table
  static std::vector<double> lookup_dist;
  
  // constructors
  Lookup() {};
  Lookup(const Data &data, const Parameters &parameters);
  
  // other mthods
  double get_data_dist(double source_lon, double source_lat, int data_index);
  
};
