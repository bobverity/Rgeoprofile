
#include <Rcpp.h>
#include <random>
#include "probability.h"
#include "misc.h"
using namespace Rcpp;

std::default_random_engine generator;
std::uniform_real_distribution<double> uniform_0_1(0.0,1.0);

//*------------------------------------------------*
// sample from uniform(a,b) distribution
double runif1(double a, double b) {
    std::uniform_real_distribution<double> uniform_a_b(a,b);
    return(uniform_a_b(generator));
}

//*------------------------------------------------*
// sample from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum) {
    double rand = pSum*uniform_0_1(generator);
    double z = 0;
    for (int i=0; i<p.size(); i++) {
        z += p[i];
        if (rand<z)
            return i+1;
    }
    return(0);
}

//*------------------------------------------------*
// sample from gamma(alpha,beta) distribution
double rgamma1(double alpha, double beta) {
    std::gamma_distribution<double> gamma_a_b(alpha,1.0/beta);
    double z = gamma_a_b(generator);
    return(z);
}

//*------------------------------------------------*
// draw from beta(alpha,beta) distribution
double rbeta1(double alpha, double beta) {
    double X1 = rgamma1(alpha,1.0);
    double X2 = rgamma1(beta,1.0);
    return(X1/(X1+X2));
}

//*------------------------------------------------*
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
    std::normal_distribution<double> normal_dist(mean,sd);
    return(normal_dist(generator));
}
