
#include <Rcpp.h>
#include "misc.h"
using namespace Rcpp;


//------------------------------------------------
// helper function for printing a single value (overloaded for int, double and string)
void print(int x) {
    Rcout << x << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}
void print(double x) {
    Rcout << x << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}
void print(std::string x) {
    Rcout << x << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}

//------------------------------------------------
// helper function for printing contents of a vector (overloaded for int, double and string)
void printVector(std::vector<int> x) {
    for (int i=0; i<int(x.size()); i++) {
        Rcout << x[i] << " ";
    }
    Rcout << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}
void printVector(std::vector<double> x) {
    for (int i=0; i<int(x.size()); i++) {
        Rcout << x[i] << " ";
    }
    Rcout << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}
void printVector(std::vector<std::string> x) {
    for (int i=0; i<int(x.size()); i++) {
        Rcout << x[i] << " ";
    }
    Rcout << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}

//------------------------------------------------
// helper function for printing contents of a matrix (overloaded for int, double and string)
void printMatrix(std::vector< std::vector<int> > M) {
    for (int i=0; i<int(M.size()); i++) {
        for (int j=0; j<int(M[i].size()); j++) {
            Rcout << M[i][j] << " ";
        }
        Rcout << "\n";
    }
    Rcout << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}
void printMatrix(std::vector< std::vector<double> > M) {
    for (int i=0; i<int(M.size()); i++) {
        for (int j=0; j<int(M[i].size()); j++) {
            Rcout << M[i][j] << " ";
        }
        Rcout << "\n";
    }
    Rcout << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}
void printMatrix(std::vector< std::vector<std::string> > M) {
    for (int i=0; i<int(M.size()); i++) {
        for (int j=0; j<int(M[i].size()); j++) {
            Rcout << M[i][j] << " ";
        }
        Rcout << "\n";
    }
    Rcout << "\n";
    R_FlushConsole();
    R_ProcessEvents();
}

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB) {
    double output = (logA<logB) ? logB + log(1+exp(logA-logB)) : logA + log(1+exp(logB-logA));
    return(output);
}

//------------------------------------------------
// take the mean of a vector
double mean(std::vector<double> &x) {
    double m = 0;
    for (int i=0; i<int(x.size()); i++) {
        m += x[i];
    }
    return(m/double(x.size()));
}

//------------------------------------------------
// round a number to defined decimal places
double round_decimal(double x, int places) {
    double scale = pow(10,places);
    double output = (x>0) ? floor(x*scale+0.5) : ceil(x*scale-0.5);
    return(output/scale);
}

//------------------------------------------------
// round a number to defined significant figures
double round_sigfig(double x, int sigfig) {
    double scale = 1;
    while (x>10) {
        x *= 0.1;
        scale *= 10;
    }
    while (x<1) {
        x *= 10;
        scale *= 0.1;
    }
    double output = scale*round_decimal(x,sigfig-1);
    return(output);
}
