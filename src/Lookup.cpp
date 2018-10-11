
#include "Lookup.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// declare static member variables

double Lookup::precision_lon;
double Lookup::precision_lat;
std::vector<double> Lookup::lookup_dist;

int Lookup::index_limit_lon_min;
int Lookup::index_limit_lon_max;
int Lookup::index_limit_lat_min;
int Lookup::index_limit_lat_max;
int Lookup::n_lon;
int Lookup::n_lat;
int Lookup::n_data;

//------------------------------------------------
// constructor for Lookup class
Lookup::Lookup(const Data &data, const Parameters &parameters) {
  
  // make copy of some parameters
  precision_lon = parameters.precision_lon;
  precision_lat = parameters.precision_lat;
  n_data = data.n;
  
  // get index limits. These are used in converting continuous lon/lat values to
  // an integer index starting at zero
  index_limit_lon_min = parameters.min_lon / precision_lon;
  index_limit_lon_max = parameters.max_lon / precision_lon;
  index_limit_lat_min = parameters.min_lat / precision_lat;
  index_limit_lat_max = parameters.max_lat / precision_lat;
  
  // populate lookup table
  n_lon = index_limit_lon_max - index_limit_lon_min + 1;
  n_lat = index_limit_lat_max - index_limit_lat_min + 1;
  
  // print size of lookup table
  //print("lookup table", sizeof(double)*n_lon*n_lat*n_data/1e6, "megabytes large");
  
  lookup_dist = vector<double>(n_lon * n_lat * n_data);
  int j = 0;
  for (int i_lon=0; i_lon<n_lon; ++i_lon) {
    for (int i_lat=0; i_lat<n_lat; ++i_lat) {
      for (int i=0; i<n_data; ++i) {
        
        // store great circle distance between this lon/lat point and the data
        double lon = (index_limit_lon_min + i_lon) * precision_lon;
        double lat = (index_limit_lat_min + i_lat) * precision_lat;
        lookup_dist[j++] = gc_dist(lon, lat, data.sentinel_lon[i], data.sentinel_lat[i]);
      }
    }
  }
  
}

//------------------------------------------------
// get distance in km of a data point (indexed by data_index) from a given 
// source lon/lat
double Lookup::get_data_dist(double source_lon, double source_lat, int data_index) {
  
  // convert lon/lat to index using precision
  int lon_index = round(source_lon / precision_lon) - index_limit_lon_min;
  int lat_index = round(source_lat / precision_lat) - index_limit_lat_min;
  
  // lookup value
  return lookup_dist[(lon_index*n_lat + lat_index)*n_data + data_index];
}