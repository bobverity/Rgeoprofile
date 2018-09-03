
#include "misc.h"
#include "hungarian.h"

using namespace std;

//------------------------------------------------
// define very large/small numbers for catching overflow/underflow problems
// OVERFLO
// UNDERFLO
// DEFINED IN HEADER

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types)
// sum
// DEFINED IN HEADER

//------------------------------------------------
// mean of vector (templated for different data types)
// mean
// DEFINED IN HEADER

//------------------------------------------------
// min of vector (templated for different data types)
// min
// DEFINED IN HEADER

//------------------------------------------------
// max of vector (templated for different data types)
// max
// DEFINED IN HEADER

//------------------------------------------------
// push back multiple values to vector
// push_back_multiple
// DEFINED IN HEADER

//------------------------------------------------
// test whether value can be found in vector
// is_in_vector
// DEFINED IN HEADER

//------------------------------------------------
// test whether two vectors have all matching values
// vectors_identical
// DEFINED IN HEADER

//------------------------------------------------
// return unique values in a contiguous sequence of integers. v_max gives the
// maximum possible value in v, which might be larger than the maximum actual
// value in v.
// Example1: v = {1 1 2 2 0 0} and v_max = 2
// Result: unique_int(v) = {1 2 0}
// Example2: v = {1 1 2 2 0 0} and v_max = 3
// Result: unique_int(v) = {1 2 0 3}
// Example3: v = {1 1 3 3 0 0} and v_max = 3
// Result: unique_int(v) = {1 3 0 2}
vector<int> unique_int(const vector<int> &v, int v_max) {
  vector<int> mask(v_max);
  vector<int> ret(v_max);
  int j = 0;
  for (int i=0; i<int(v.size()); i++) {
    if (mask[v[i]]==0) {
      ret[j] = v[i];
      mask[v[i]] = 1;
      j++;
    }
  }
  for (int k=0; k<v_max; k++) {
    if (mask[k]==0) {
      ret[j] = k;
      j++;
    }
  }
  return ret;
}

//------------------------------------------------
// return order of unique values in a contiguous sequence of integers. v_max 
// gives the maximum possible value in v, which might be larger than the maximum
// actual value in v.
// Example1: v = {1 1 2 2 0 0} and v_max = 2
// Result: order_unique_int(v) = {2 0 1}
// Example2: v = {1 1 2 2 0 0} and v_max = 3
// Result: order_unique_int(v) = {2 0 1 3}
// Example3: v = {1 1 3 3 0 0} and v_max = 3
// Result: order_unique_int(v) = {2 0 3 1}
vector<int> order_unique_int(const vector<int> &v, int v_max) {
  vector<int> ret(v_max, -1);
  int j = 0;
  for (int i=0; i<int(v.size()); i++) {
    if (ret[v[i]]<0) {
      ret[v[i]] = j;
      j++;
    }
  }
  for (int k=0; k<v_max; k++) {
    if (ret[k]<0) {
      ret[k] = j;
      j++;
    }
  }
  return ret;
}

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double log_sum(double logA, double logB) {
  if (logA-logB > 100) {
    return logA;
  } else if (logB-logA > 100) {
    return logB;
  }
  double ret = (logA<logB) ? logB + log(1+exp(logA-logB)) : logA + log(1+exp(logB-logA));
  return ret;
}

//------------------------------------------------
// call Hungarian algorithm for binding best matching in a linear sum assigment problem
// [[Rcpp::export]]
Rcpp::List call_hungarian_cpp(Rcpp::List args) {
  
  // objects for calling Hungarian algorithm
  vector<vector<double>> cost_mat = rcpp_to_mat_double(args["cost_mat"]);
  int n = cost_mat.size();
  vector<int> edges_left(n);
  vector<int> edges_right(n);
  vector<int> blocked_left(n);
  vector<int> blocked_right(n);
  vector<int> best_perm = hungarian(cost_mat, edges_left, edges_right, blocked_left, blocked_right);
  
  // return
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("best_matching") = best_perm);
  return ret;
}

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
// print
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
// print_vector
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
// print_matrix
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
// print_array
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of an Rcpp numeric vector
void rcpp_print_vector(Rcpp::NumericVector &x) {
  for (int i=0; i<x.length(); i++) {
    Rcpp::Rcout << x[i] << " ";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}
void rcpp_print_vector(Rcpp::IntegerVector &x) {
  for (int i=0; i<x.length(); i++) {
    Rcpp::Rcout << x[i] << " ";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of an Rcpp numeric matrix
void rcpp_print_matrix(Rcpp::NumericMatrix &x) {
  for (int i=0; i<x.nrow(); i++) {
    for (int j=0; j<x.ncol(); j++) {
      Rcpp::Rcout << x(i,j) << " ";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}
void rcpp_print_matrix(Rcpp::IntegerMatrix &x) {
  for (int i=0; i<x.nrow(); i++) {
    for (int j=0; j<x.ncol(); j++) {
      Rcpp::Rcout << x(i,j) << " ";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(string title, int n) {
  Rcpp::Rcout << title << " ";
  for (int i=0; i<n; i++) {
    Rcpp::Rcout << "*";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// print "foo", with option number e.g. "foo2"
void foo(int n) {
  if (n==0) {
    Rcpp::Rcout << "foo\n";
  } else {
    Rcpp::Rcout << "foo" << n << "\n";
  }
  R_FlushConsole();
}

//------------------------------------------------
// print "bar", with option number e.g. "bar2"
void bar(int n) {
  if (n==0) {
    Rcpp::Rcout << "bar\n";
  } else {
    Rcpp::Rcout << "bar" << n << "\n";
  }
  R_FlushConsole();
}

//------------------------------------------------
// print "foobar", with option number e.g. "foobar2"
void foobar(int n) {
  if (n==0) {
    Rcpp::Rcout << "foobar\n";
  } else {
    Rcpp::Rcout << "foobar" << n << "\n";
  }
  R_FlushConsole();
}

//------------------------------------------------
// analogue of R function seq() for integers
vector<int> seq_int(int from, int to, int by) {
  if (by==0) {
    Rcpp::stop("called seq_int with by=0");
  }
  int n = floor((to-from)/double(by)) + 1;
  vector<int> ret(n,from);
  for (int i=1; i<n; i++) {
    from += by;
    ret[i] = from;
  }
  return ret;
}

//------------------------------------------------
// converts input from Rcpp::List format to bool format.
int rcpp_to_bool(SEXP x) {
  return Rcpp::as<bool>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to int format.
int rcpp_to_int(SEXP x) {
  return Rcpp::as<int>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to double format.
double rcpp_to_double(SEXP x) {
  return Rcpp::as<double>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to string format.
string rcpp_to_string(SEXP x) {
  return Rcpp::as<string>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<bool> format.
vector<bool> rcpp_to_vector_bool(SEXP x) {
  return Rcpp::as<vector<bool> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
vector<int> rcpp_to_vector_int(SEXP x) {
  return Rcpp::as<vector<int> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
vector<double> rcpp_to_vector_double(SEXP x) {
  return Rcpp::as<vector<double> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<string> format.
vector<string> rcpp_to_vector_string(SEXP x) {
  return Rcpp::as<vector<string> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<bool>> format.
vector< vector<bool> > rcpp_to_mat_bool(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<bool> > x_mat(nrow);
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<bool> >(x[i]);
  }
  return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
vector< vector<int> > rcpp_to_mat_int(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<int> > x_mat(nrow);
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<int> >(x[i]);
  }
  return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
vector< vector<double> > rcpp_to_mat_double(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<double> > x_mat(nrow);
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<double> >(x[i]);
  }
  return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
vector< vector< vector<int> > > rcpp_to_array_int(Rcpp::List x) {
  int n1 = int(x.size());
  vector< vector< vector<int> > > ret(n1);
  for (int i=0; i<n1; i++) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector< vector<int> >(n2);
    for (int j=0; j<n2; j++) {
      ret[i][j] = Rcpp::as<vector<int> >(x_i[j]);
    }
  }
  return ret;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
vector< vector< vector<double> > > rcpp_to_array_double(Rcpp::List x) {
  int n1 = int(x.size());
  vector< vector< vector<double> > > ret(n1);
  for (int i=0; i<n1; i++) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector< vector<double> >(n2);
    for (int j=0; j<n2; j++) {
      ret[i][j] = Rcpp::as<vector<double> >(x_i[j]);
    }
  }
  return ret;
}

// --------------------------------
// Bhaskara I's approximation to sine, cosine, and tanjent functions. Sine and
// cosine approximations have a maximum deviation from the true value of approx.
// 0.00163
double quick_cos(double x) {
  
  // bring x into the range [0,2*PI]
  while (x<0) {
    x += 2*PI;
  }
  while (x>2*PI) {
    x -= 2*PI;
  }
  
  // divide [0,2*PI] into 4 intervals. Use approximation on each interval
  if (x<PI/2) {
    return (PI_sq-4*x*x)/(PI_sq+x*x);
  } else if (x<PI) {
    x = PI - x;
    return -(PI_sq-4*x*x)/(PI_sq+x*x);
  } else if (x<(3*PI/2)) {
    x -= PI;
    return -(PI_sq-4*x*x)/(PI_sq+x*x);
  } else {
    x = 2*PI - x;
    return (PI_sq-4*x*x)/(PI_sq+x*x);
  }
}
double quick_sin(double x) {
  return quick_cos(x-PI/2);
}
double quick_tan(double x) {
  return quick_sin(x)/quick_cos(x);
}

// --------------------------------
// get great circle distance between two points in lon/lat coordinates
double gc_dist(double lon0, double lat0, double lon1, double lat1) {
  
  // check for exact equality of points
  if (lon0 == lon1 && lat0 == lat1) {
    return 0;
  }
  
  // convert input arguments to radians
  lon0 *= 2*PI/360;
  lat0 *= 2*PI/360;
  lon1 *= 2*PI/360;
  lat1 *= 2*PI/360;
  
  // calculate great circle angle. Use temporary variable to avoid acos(>1) or
  // acos(<0), which can occur due to underflow issues
  double tmp = sin(lat0)*sin(lat1) + cos(lat0)*cos(lat1)*cos(lon1-lon0);
  tmp = (tmp > 1.0) ? 1.0 : tmp;
  tmp = (tmp < 0.0) ? 0.0 : tmp;
  double gc_angle = acos(tmp);
  
  // convert gc_angle to great circle distance using radius of earth (km)
  double gc_dist = gc_angle * EARTH_RAD_KM;
  
  return gc_dist;
}


