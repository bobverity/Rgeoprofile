
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing data
class Data {
  
public:
  
  static std::vector<double> sentinel_lon;
  static std::vector<double> sentinel_lat;
  static std::vector<int> sentinel_counts;
  static int n;
  
  // constructors
  Data() {};
  Data(const Rcpp::List &args);
  
};
