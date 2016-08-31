
#ifndef __RgeoProfile__main__
#define __RgeoProfile__main__

using namespace std;

// carry out main MCMC
Rcpp::List C_geoMCMC(Rcpp::List data, Rcpp::List params);

// update group allocation
void updateGroup(int &i, int &n, vector<double> &data_x, vector<double> &data_y, vector<int> &group, vector<int> &freqs, vector<double> &sum_x, vector<double> &sum_y, vector<double> &sumSquared_x, vector<double> &sumSquared_y, int &nextGroup, int &uniqueGroups, vector<double> &mu_postVar, vector<double> &mu_postMean_x, vector<double> &mu_postMean_y, vector<double> &mu_postDraw_x, vector<double> &mu_postDraw_y, double &sigma2, double &tau2, vector<double> &probVec, vector<double> &logProbVec, double &alpha);

// solve label switching problem
void solveLabelSwitching(int &n, vector<int> &group, vector< vector<int> > &groupMat, vector<int> &bestPerm, vector<int> &bestPermOrder, vector<int> &group_reorder, vector<int> &freqs, vector<double> &sum_x, vector<double> &sum_y, vector<double> &sumSquared_x, vector<double> &sumSquared_y, vector<int> &freqs_reorder, vector<double> &sum_x_reorder, vector<double> &sum_y_reorder, vector<double> &sumSquared_x_reorder, vector<double> &sumSquared_y_reorder, int &nextGroup);

#endif