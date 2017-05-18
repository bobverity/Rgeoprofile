
#ifndef __RgeoProfile__probability__
#define __RgeoProfile__probability__

//------------------------------------------------
// sample from uniform(a,b) distribution
double runif1(double a=0, double b=1.0);

//------------------------------------------------
// sample from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum=1);

//------------------------------------------------
// sample from gamma(alpha,beta) distribution
double rgamma1(double shape, double rate);

//------------------------------------------------
// draw from beta(alpha,beta) distribution
double rbeta1(double alpha, double beta);

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd);

#endif
