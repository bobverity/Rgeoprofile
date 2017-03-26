
#ifndef __RgeoProfile__misc__
#define __RgeoProfile__misc__


//------------------------------------------------
// helper function for printing a single value (overloaded for int, double and string)
void print(int x);
void print(double x);
void print(std::string x);

//------------------------------------------------
// helper function for printing contents of a vector (overloaded for int, double and string)
void printVector(std::vector<int> x);
void printVector(std::vector<double> x);
void printVector(std::vector<std::string> x);

//------------------------------------------------
// helper function for printing contents of a matrix (overloaded for int, double and string)
void printMatrix(std::vector< std::vector<int> > M);
void printMatrix(std::vector< std::vector<double> > M);
void printMatrix(std::vector< std::vector<std::string> > M);

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB);

//------------------------------------------------
// take the mean of a vector
double mean(std::vector<double> &x);

//------------------------------------------------
// round a number to defined decimal places
double round_decimal(double x, int places=0);

//------------------------------------------------
// round a number to defined significant figures
double round_sigfig(double x, int sigfig=1);

#endif
