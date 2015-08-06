
#ifndef __RgeoProfile__main__
#define __RgeoProfile__main__

using namespace Rcpp;

// carry out main MCMC
List C_geoMCMC(List data, List params);

// update group allocation
void updateGroup(int &i, int &n, std::vector<double> &data_x, std::vector<double> &data_y, std::vector<int> &group, std::vector<int> &freqs, std::vector<double> &sum_x, std::vector<double> &sum_y, int &nextGroup, int &uniqueGroups, std::vector<double> &postVar, std::vector<double> &postMean_x, std::vector<double> &postMean_y, double &priorMean_x, double &priorMean_y, double &sigma2, double &tau2, std::vector<double> &probVec, std::vector<double> &logProbVec, double &alpha);

// solve label switching problem
void solveLabelSwitching(int &n, std::vector<int> &group, std::vector< std::vector<int> > &groupMat, std::vector<int> &bestPerm, std::vector<int> &bestPermOrder, std::vector<int> &group_reorder, std::vector<int> &freqs, std::vector<double> &sum_x, std::vector<double> &sum_y, std::vector<int> &freqs_reorder, std::vector<double> &sum_x_reorder, std::vector<double> &sum_y_reorder, int &nextGroup);

#endif