
#pragma once

#include <random>

//------------------------------------------------
// Hungarian algorithm for solving linear-sum assignment problem (further notes in body file)
void augment_left(std::vector<int> &ret, int i, std::vector< std::vector<double> > &m, std::vector<int> &edges_right, std::vector<int> &blocked_left, std::vector<int> &blocked_right);
void augment_right(std::vector<int> &ret, int j, std::vector< std::vector<double> > &m, std::vector<int> &edges_right, std::vector<int> &blocked_left, std::vector<int> &blocked_right);
std::vector<int> hungarian(std::vector< std::vector<double> > &m, std::vector<int> &edges_left, std::vector<int> &edges_right, std::vector<int> &blocked_left, std::vector<int> &blocked_right, const int max_reps=int(1e6));

