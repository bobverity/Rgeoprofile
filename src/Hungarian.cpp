
#include <Rcpp.h>
#include <random>
#include "Hungarian.h"
#include "misc.h"
using namespace Rcpp;

//*------------------------------------------------*
// transpose matrix
std::vector< std::vector<double> > transpose(std::vector< std::vector<double> > &M) {
    std::vector< std::vector<double> > tM(M[0].size());
    for (int i=0; i<M[0].size(); i++) {
        tM[i] = std::vector<double>(M.size());
        for (int j=0; j<M.size(); j++) {
            tM[i][j] = M[j][i];
        }
    }
    return(tM);
}

//*------------------------------------------------*
// for each row in a matrix, subtract smallest element in that row
void subtractRows(std::vector< std::vector<double> > &M) {
    double minVal;
    for (int i=0; i<M.size(); i++) {
        minVal = *min_element(begin(M[i]),end(M[i]));
        for (unsigned int j=0; j<M[i].size(); j++) {
            M[i][j] -= minVal;
        }
    }
}

//*------------------------------------------------*
// for each column in a matrix, subtract smallest element in that column
void subtractCols(std::vector< std::vector<double> > &M) {
    std::vector< std::vector<double> > tM = transpose(M);
    double minVal;
    for (int i=0; i<tM.size(); i++) {
        minVal = *std::min_element(std::begin(tM[i]),std::end(tM[i]));
        for (int j=0; j<tM[i].size(); j++) {
            M[j][i] -= minVal;
        }
    }
}

//*------------------------------------------------*
// the functions augmentLeft and augmentRight work together to find an augmented path. They call each other, which normally could lead to an infinite recursion, but this is avoided as eventually either an augmented path will be found or no more moves will be possible. Return full path, or -1 if no path found.
std::vector<int> augmentLeft(int i, std::vector< std::vector<double> > &M, std::vector<int> &edgesRight, std::vector<int> &blockedLeft, std::vector<int> &blockedRight) {
    
    blockedLeft[i] = 1;
    std::vector<int> output(1,-1);
    
    // search all unmatched edges
    for (int j=0; j<M.size(); j++) {
        if (M[i][j]==0 && edgesRight[j]!=i && blockedRight[j]==0) {
            
            // if edge leads to augmented path then add current node to path and return
            output = augmentRight(j, M, edgesRight, blockedLeft, blockedRight);
            if (output[0]>=0) {
                output.push_back(i);
                return(output);
            }
        }
    }
    
    // if no more moves then return -1
    return(output);
}

std::vector<int> augmentRight(int j, std::vector< std::vector<double> > &M, std::vector<int> &edgesRight, std::vector<int> &blockedLeft, std::vector<int> &blockedRight) {
    
    blockedRight[j] = 1;
    std::vector<int> output(1);
    
    // if node j is unmatched then return j as start of augmented path
    if (edgesRight[j]<0) {
        output[0] = j;
        return(output);
    }
    
    // otherwise continue chain of augmenting
    output = augmentLeft(edgesRight[j], M, edgesRight, blockedLeft, blockedRight);
    if (output[0]>=0) {
        output.push_back(j);
    }
    return(output);
}

//*------------------------------------------------*
// carry out Hungarian algorithm to find best matching given cost matrix M
std::vector<int> hungarian(std::vector< std::vector<double> > &M) {
    
    int maxReps= 1000000;
    int n = int(M.size());
    
    // initialise assignment objects
    std::vector<int> edgesLeft;
    std::vector<int> edgesRight;
    int numberAssigned;
    
    // search for solution until maxReps reached
    for (int rep=0; rep<maxReps; rep++) {
        
        // zero assignment objects
        edgesLeft = std::vector<int>(n,-1);
        edgesRight = std::vector<int>(n,-1);
        numberAssigned = 0;
        
        // subtract smallest element from all rows and columns
        subtractRows(M);
        subtractCols(M);
        
        // generate an initial matching
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (M[i][j]==0 && edgesRight[j]<0) {
                    edgesLeft[i] = j;
                    edgesRight[j] = i;
                    numberAssigned++;
                    break;
                }
            }
        }
        
        // if this matching is perfect then we are done
        if (numberAssigned==n)
            return(edgesLeft);
        
        // continue augmenting paths until no more possible
        std::vector<int> blockedLeft(n,0);
        std::vector<int> blockedRight(n,0);
        bool continueAugmenting = true;
        while (continueAugmenting==true) {
            continueAugmenting = false;
            
            // search all unmatched nodes
            for (int i=0; i<n; i++) {
                if (edgesLeft[i]<0) {
                    
                    // attempt to find augmented path
                    blockedLeft = std::vector<int>(n,0);
                    blockedRight = std::vector<int>(n,0);
                    std::vector<int> path = augmentLeft(i, M, edgesRight, blockedLeft, blockedRight);
                    
                    // if successful then augment
                    if (path[0]>=0) {
                        continueAugmenting = true;
                        numberAssigned ++;
                        for (int j=0; j<(path.size()/2); j++) {
                            edgesLeft[path[j*2+1]] = path[j*2];
                            edgesRight[path[j*2]] = path[j*2+1];
                        }
                        
                        // if best matching found then finish
                        if (numberAssigned==n) {
                            return(edgesLeft);
                        }
                    }
                }
            }
        }
        
        // find minimum value in cost matrix, looking at all elements in which neither the row or the column is part of the minimum vertex cover
        double minVal = -log(double(0));
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (blockedLeft[i]==1 && blockedRight[j]==0 && M[i][j]<minVal) {
                    minVal = M[i][j];
                }
            }
        }
        
        // add or subtract this value to cost matrix as required
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (blockedLeft[i]==1 && blockedRight[j]==0) {
                    M[i][j] -= minVal;
                }
                if (blockedLeft[i]==0 && blockedRight[j]==1) {
                    M[i][j] += minVal;
                }
            }
        }
        
        // at this point we have a new cost matrix and can repeat the process from the top
        
    } // rep loop
    
    Rcout << "Error in Hungarian algorithm\n";
    exit(1);
    
}


