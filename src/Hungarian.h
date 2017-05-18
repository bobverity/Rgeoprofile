
#ifndef __RgeoProfile__Hungarian__
#define __RgeoProfile__Hungarian__

// transpose matrix
std::vector< std::vector<double> > transpose(std::vector< std::vector<double> > &M);

// for each row in a matrix, subtract smallest element in that row
void subtractRows(std::vector< std::vector<double> > &M);

// for each column in a matrix, subtract smallest element in that column
void subtractCols(std::vector< std::vector<double> > &M);

// the functions augmentLeft and augmentRight work together to find an augmented path. They call each other, which normally could lead to an infinite recursion, however, this is avoided as eventually an augmented path will be found, or no more moves will be possible. Return full path, or -1 if no path found.
std::vector<int> augmentLeft(int i, std::vector< std::vector<double> > &M, std::vector<int> &edgesRight, std::vector<int> &blockedLeft, std::vector<int> &blockedRight);
std::vector<int> augmentRight(int j, std::vector< std::vector<double> > &M, std::vector<int> &edgesRight, std::vector<int> &blockedLeft, std::vector<int> &blockedRight);

// carry out Hungarian algorithm to find best matching given cost matrix M
std::vector<int> hungarian(std::vector< std::vector<double> > &M);

#endif
