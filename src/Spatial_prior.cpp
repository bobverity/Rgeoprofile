
#include "Spatial_prior.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// declare static member variables

vector<bool> Spatial_prior::spatial_prior_mask;

//------------------------------------------------
// constructor for Parameters class
Spatial_prior::Spatial_prior(const Rcpp::List &args) {
  
  // get number of cells in each dimension
  vector<double> spatial_prior_values = rcpp_to_vector_double(args["spatial_prior_values"]);
  
  // populate mask
  spatial_prior_mask = vector<bool>(n_lat*n_lon, false);
  for (int i=0; i<(n_lat*n_lon); ++i) {
    spatial_prior_mask[i] = (spatial_prior_values[i] != 0) ? true : false;
  }
  
}

//------------------------------------------------
// get value
double Spatial_prior::get_value(double lon, double lat) {
  
  // convert lon/lat to index using precision
  int lon_index = floor(lon / res_lon) - index_lon_min;
  int lat_index = floor(lat / res_lat) - index_lat_min;
  
  // lookup value
  return spatial_prior_mask[(n_lat-lat_index)*n_lon + lon_index];
  //return spatial_prior_mask[lon_index*n_lat + lat_index];
  
}