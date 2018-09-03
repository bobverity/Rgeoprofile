
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// define pi squared
#define PI_sq 9.86960440
#define LOG_PI 1.14472989

//------------------------------------------------
// define very large/small numbers for catching overflow/underflow problems
#define OVERFLO   1e100
#define UNDERFLO   1e-100

//------------------------------------------------
// define radius of earth in km
#define EARTH_RAD_KM   6371

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types).
template<class TYPE>
TYPE sum(std::vector<TYPE> &x) {
    TYPE output = 0;
    for (int i=0; i<int(x.size()); i++) {
        output += x[i];
    }
    return output;
}

//------------------------------------------------
// mean of vector (templated for different data types)
template<class TYPE>
double mean(std::vector<TYPE> &x) {
    return sum(x)/double(x.size());
}

//------------------------------------------------
// min of vector (templated for different data types)
template<class TYPE>
TYPE min(std::vector<TYPE> x) {
    return *min_element(x.begin(), x.end());
}

//------------------------------------------------
// max of vector (templated for different data types)
template<class TYPE>
TYPE max(std::vector<TYPE> x) {
    return *max_element(x.begin(), x.end());
}

//------------------------------------------------
// push back multiple values to vector
template<class TYPE>
void push_back_multiple(std::vector<TYPE> &lhs, std::vector<TYPE> &rhs) {
    lhs.insert(lhs.end(), rhs.begin(), rhs.end());
}

//------------------------------------------------
// test whether value can be found in vector
template<class TYPE>
bool is_in_vector(TYPE x, std::vector<TYPE> &v) {
  bool ret = false;
  for (int i=0; i<int(v.size()); i++) {
    if (v[i] == x) {
      ret = true;
      break;
    }
  }
  return ret;
}

//------------------------------------------------
// test whether two vectors have all matching values
template<class TYPE>
bool vectors_identical(std::vector<TYPE> &v1, std::vector<TYPE> &v2) {
  for (int i=0; i<int(v1.size()); i++) {
    if (v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

//------------------------------------------------
// return unique values in a contiguous sequence of integers. v_max gives the
// maximum possible value in v, which might be larger than the maximum actual
// value in v.
std::vector<int> unique_int(const std::vector<int> &v, int v_max);

//------------------------------------------------
// return order of unique values in a contiguous sequence of integers. v_max 
// gives the maximum possible value in v, which might be larger than the maximum
// actual value in v.
std::vector<int> order_unique_int(const std::vector<int> &v, int v_max);

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double log_sum(double logA, double logB);

//------------------------------------------------
// call Hungarian algorithm for binding best matching in a linear sum assigment problem
Rcpp::List call_hungarian_cpp(Rcpp::List args);

//------------------------------------------------
// helper function for printing a single value or series of values (templated for different data types)
template<class TYPE>
void print(TYPE x) {
    Rcpp::Rcout << x << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2>
void print(TYPE1 x1, TYPE2 x2) {
    Rcpp::Rcout << x1 << " " << x2 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4, TYPE5 x5) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5, class TYPE6>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4, TYPE5 x5, TYPE6 x6) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
template<class TYPE>
void print_vector(std::vector<TYPE> &x) {
    for (int i=0; i<x.size(); i++) {
        Rcpp::Rcout << x[i] << " ";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
template<class TYPE>
void print_matrix(std::vector< std::vector<TYPE> > &x) {
    for (int i=0; i<x.size(); i++) {
        for (int j=0; j<x[i].size(); j++) {
            Rcpp::Rcout << x[i][j] << " ";
        }
        Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
template<class TYPE>
void print_array(std::vector< std::vector< std::vector<TYPE> > > &x) {
    for (int i=0; i<x.size(); i++) {
        Rcpp::Rcout << "--- slice " << i+1 << " ---\n";
        for (int j=0; j<x[i].size(); j++) {
            for (int k=0; k<x[i][j].size(); k++) {
                Rcpp::Rcout << x[i][j][k] << " ";
            }
            Rcpp::Rcout << "\n";
        }
        Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of an Rcpp numeric vector
void rcpp_print_vector(Rcpp::NumericVector &x);
void rcpp_print_vector(Rcpp::IntegerVector &x);

//------------------------------------------------
// helper function for printing contents of an Rcpp numeric matrix
void rcpp_print_matrix(Rcpp::NumericMatrix &x);
void rcpp_print_matrix(Rcpp::IntegerMatrix &x);

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(std::string title="", int n=10);

//------------------------------------------------
// print "foo", with option number e.g. "foo2"
void foo(int n=0);

//------------------------------------------------
// print "bar", with option number e.g. "bar2"
void bar(int n=0);

//------------------------------------------------
// print "foobar", with option number e.g. "foobar2"
void foobar(int n=0);

//------------------------------------------------
// analogue of R function seq() for integers
std::vector<int> seq_int(int from, int to, int by=1);

//------------------------------------------------
// converts input from Rcpp::List format to bool format.
int rcpp_to_bool(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to int format.
int rcpp_to_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to double format.
double rcpp_to_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to string format.
std::string rcpp_to_string(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<bool> format.
std::vector<bool> rcpp_to_vector_bool(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
std::vector<int> rcpp_to_vector_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
std::vector<double> rcpp_to_vector_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<string> format.
std::vector<std::string> rcpp_to_vector_string(SEXP x);

// converts input from Rcpp::List format to vector<vector<bool>> format.
std::vector< std::vector<bool> > rcpp_to_mat_bool(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
std::vector< std::vector<int> > rcpp_to_mat_int(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
std::vector< std::vector<double> > rcpp_to_mat_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
std::vector< std::vector< std::vector<double> > > rcpp_to_array_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
std::vector< std::vector< std::vector<int> > > rcpp_to_array_int(Rcpp::List x);

// --------------------------------
// Bhaskara I's approximation to sine, cosine, and tanjent functions
double quick_cos(double x);
double quick_sin(double x);
double quick_tan(double x);

// --------------------------------
// get great circle distance between two points in lon/lat coordinates
double gc_dist(double lon0, double lat0, double lon1, double lat1);



