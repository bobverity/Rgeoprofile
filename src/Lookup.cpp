
#include "Lookup.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// declare static member variables

vector<double> Lookup::lookup_dist;

//------------------------------------------------
// constructor for Lookup class
Lookup::Lookup() {
  
  // print size of lookup table
  //print("lookup table", sizeof(double)*n_lon*n_lat*n/1e6, "megabytes large");
  
  // populate lookup table
  lookup_dist = vector<double>(n_lon * n_lat * n);
  int j = 0;
  for (int i_lon=0; i_lon<n_lon; ++i_lon) {
    for (int i_lat=0; i_lat<n_lat; ++i_lat) {
      for (int i=0; i<n; ++i) {
        
        // store great circle distance between this lon/lat point and the data
        double lon = (index_lon_min + i_lon + 0.5)*res_lon;
        double lat = (index_lat_min + i_lat + 0.5)*res_lat;
        lookup_dist[j++] = gc_dist(lon, lat, sentinel_lon[i], sentinel_lat[i]);
      }
    }
  }
  
}

//------------------------------------------------
// get distance in km of a data point (indexed by data_index) from a given 
// source lon/lat
double Lookup::get_data_dist(double source_lon, double source_lat, int data_index) {
  
  // convert lon/lat to index
  int lon_index;
  int lat_index;
  if (source_lon > 0) {
    lon_index = floor(source_lon / res_lon) - index_lon_min;
  } else {
    lon_index = ceil(source_lon / res_lon) - index_lon_min;
  }
  if (source_lat > 0) {
    lat_index = floor(source_lat / res_lat) - index_lat_min;
  } else {
    lat_index = ceil(source_lat / res_lat) - index_lat_min;
  }
  
  // lookup value
  return lookup_dist[(lon_index*n_lat + lat_index)*n + data_index];
}