
#include "hungarian.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// HUNGARIAN ALGORITHM
// For solving linear-sum assignment problem to return lowest-cost assignment.
//
// ## Motivation
// Imagine a set of n nodes on the left representing workers, and n nodes on the
// right representing tasks to be done. Edges between nodes on the left and 
// nodes on the right indicate that a task can be done by a worker, and the 
// weight of the edge represents the cost incurred in doing this task. We want 
// to find the assignment such that every worker is assigned one task, and such
// that the overall cost is minimized.
//
//
// ## Simple case
// We can visualise the cost of all edges in a cost matrix (denoted m), with 
// rows representing workers and columns representing tasks. Take the following 
// 3-by-3 matrix
// 
// 4 2 5
// 1 4 7
// 6 6 3
// 
// When comparing the edges of a node, all that matters is the relative cost. 
// Therefore we can subtract the smallest value in every row and column of m to 
// arrive at an equivalent cost matrix, where best matchings are now represented
// by 0s. The new matrix becomes:
// 
// 2 0 3
// 0 3 6
// 3 3 0
// 
// We now create two vectors, edges_left and edges_right, both of length 3 (the 
// dimension of m). edges_left will store which task each worker is assigned to,
// and edges_right will store which worker each task is assigned to. Both
// vectors are initialised with -1 values to indicate that no assigmnet has been
// made.
// We loop through all rows and columns of m. If m[i][j] is equal to 0, and if
// edges_right[j] is currently unassigned then we create a new edge between node
// i on the left and node j on the right as follows:
//
// edges_left[i] = j
// edges_right[j] = i
//
// In the simple case above, this create the following perfect matching:
//
// edges_left = {1, 0, 2}
// edges_right = {1, 0, 2}
//
// (note that these vectors will not always contain the same values). The vector
// edges_left is returned as the solution.
//
// 
// ## Augmenting paths
// Now imagine that m, after subtracting the minimum of rows and columns, looks
// like this:
// 
// 0 9 9
// 9 0 0
// 9 0 9
//
// The first pass will produce the following imperfect matching:
//
// edges_left = {0, 1, -1}
// edges_right = {0, 1, -1}
//
// We can draw this matching (rather badly) as follows, where double lines
// indicate the chosen matching, and single lines indicate alternative possible
// zero-cost matchings.
//
// * ========== *
// * ========== *
//   \___  ___/
//       \/
//   ____/\____
// *            *
//
// Note that we could improve the matching if the second worker was assigned the
// third task and vice versa. In general, we can imagine starting from an 
// unmatched worker with a zero-cost edge (worker 3), traverse this edge to a 
// task with a matched edge (task 2), traverse this edge to worker with an 
// unmatched zero-cost edge (worker 1) and traverse this edge to a task (task 
// 3), which as no matched worker, therefore we stop and terminate the path. 
// This path alternates between unmatched and matched edges, and always starts 
// and ends on unmatched edges, meaning there is one more unmatched than matched
// edge. We can therefore swap matched and unmatched edges along this path, and 
// we will always augment our number of matches by 1. In this example, a single
// augmenting path brings us to the complete solution, although in general
// multiple augmenting paths may be needed, and even then it is possible that a
// solution will not be found.
//
//
// ## Adjusting the cost matrix
// If there are no more augmenting paths, and a complete assignment has still
// not been found, then we need to adjust the cost matrix and try again. In
// simple terms, we want to lock-in the zero-cost edges that we have found, and
// then find the second minimum cost among the remaining choices.
// Take the following cost matrix:
// 
// 0 2 3
// 0 1 4
// 0 0 0
//
// Here the initial matching will produce
//
// edges_left = {0, -1, 1}
// edges_right = {0, 2, -1}
//
// Augmenting paths will not improve upon this solution. Therefore, we want to 
// lock-in the current zero-cost edges in the most parsimonious way, i.e. 
// blocking off as few nodes as possible. Thinking in terms of the cost matrix, 
// this is equivalent to drawing lines through rows and columns that cross out 
// all 0s, while drawing as few lines as possible. Clearly in this case the best
// option will be to cross of the first column and the final row. More 
// generally, which rows and columns to cross off can be determined from the 
// previous augmented paths step. We simply mark every node that was visited 
// when searching for augmented paths - in this example this will be the second 
// and first nodes on the left, and the first node on the right. We then cross 
// out all marked nodes on the right, and all *unmarked* nodes on the left. This
// is equivalent to crossing out the first column and the final row of the cost 
// matrix. We then calculate the minimum value of all remaining (un-crossed-out)
// elements. We subtract this value from the unmarked elements, and add this 
// value to any elements that are doubly crossed out (in this example the
// bottom-left zero). This yeilds the following cost matrix:
// 
// 0 1 2
// 0 0 3
// 1 0 0
//
// We then repeat the entire procedure of searching for a matching and
// augmenting paths etc. on this new cost matrix, continuing until we have found
// a complete matching. The matching returned by this procedure is guaranteed to
// be the lowest-cost complete matching.
//
//
// ## Implementation details
// The main part of the algorithm, in hungarian(), can be easily understaood 
// from the description above. The only notable difference is that the simple 
// identity matching (i.e. a cost matrix with diagonal zeros) is searched for 
// manually first, as in some applications this matching will be extremely 
// common and so a manual search is computationally efficient. The most complex 
// part of the implementation are the functions augment_left() and 
// augment_right(). These functions call each other recursively, stepping 
// forward through the augmenting path. If the path fails to terminate on an 
// unmatched edge with zero weight, then the function returns a value -1 in the 
// ret vector. This -1 value is then propagated backwards through the path back 
// to the start, indicating that no augmented path was found. An infinite 
// recursion in these functions is avoided, as eventually either an augmented 
// path will be found or the recursion will terminate. Note that hungarian()
// contains a cut-out limit of max_reps - if no complete matching has been found
// within this many repeats of updating the cost matrix then the function stops
// using Rcpp::stop().

//------------------------------------------------

void augment_left(vector<int> &ret, int i, vector< vector<double> > &m, vector<int> &edges_right, vector<int> &blocked_left, vector<int> &blocked_right) {
  
  // block node i on the left
  blocked_left[i] = 1;
  
  // search all unmatched edges
  for (int j=0; j<int(m.size()); j++) {
    if (m[i][j]==0 && blocked_right[j]==0) {
      
      // if edge leads to augmented path then add current node to path and return
      augment_right(ret, j, m, edges_right, blocked_left, blocked_right);
      if (ret[0]>=0) {
        ret.push_back(i);
        return;
      }
    }
  }
  
  // if no more moves then return as is
  return;
}

//------------------------------------------------

void augment_right(vector<int> &ret, int j, vector< vector<double> > &m, vector<int> &edges_right, vector<int> &blocked_left, vector<int> &blocked_right) {
  
  // block node j on the right
  blocked_right[j] = 1;
  
  // if node j is unmatched then return j as start of augmented path
  if (edges_right[j]<0) {
    ret.back() = j;
    return;
  }
  
  // otherwise continue chain of augmenting
  augment_left(ret, edges_right[j], m, edges_right, blocked_left, blocked_right);
  if (ret[0]>=0) {
    ret.push_back(j);
  }
  
}

//------------------------------------------------

vector<int> hungarian(vector< vector<double> > &m, vector<int> &edges_left, vector<int> &edges_right, vector<int> &blocked_left, vector<int> &blocked_right, const int max_reps) {
  
  // initialise objects
  int n = m.size();
  vector<double> min_col(n);
  
  // search for solution until max_reps reached
  for (int rep=0; rep<max_reps; rep++) {
    
    // reset search objects
    fill(edges_left.begin(), edges_left.end(), -1);
    fill(edges_right.begin(), edges_right.end(), -1);
    int number_assigned = 0;
    
    // subtract smallest element from all rows and columns
    min_col = m[0];
    for (int i=0; i<n; i++) {
      double min_row = min(m[i]);
      for (int j=0; j<n; j++) {
        m[i][j] -= min_row;
        if (m[i][j]<min_col[j]) {
          min_col[j] = m[i][j];
        }
      }
    }
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        m[i][j] -= min_col[j];
      }
    }
    
    // check for simple identity matching (most common case)
    bool identity_matching = true;
    for (int i=0; i<n; i++) {
      if (m[i][i]!=0) {
        identity_matching = false;
        break;
      }
    }
    if (identity_matching) {
      return seq_int(0,n-1);
    }
    
    // generate an initial matching
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (m[i][j]==0 && edges_right[j]<0) {
          edges_left[i] = j;
          edges_right[j] = i;
          number_assigned ++;
          break;
        }
      }
    }
    
    // if this matching is perfect then we are done
    if (number_assigned==n) {
      return edges_left;
    }
    
    // continue augmenting paths until no more possible
    bool continue_augmenting = true;
    while (continue_augmenting) {
      continue_augmenting = false;
      
      // zero blocking vectors
      fill(blocked_left.begin(), blocked_left.end(), 0);
      fill(blocked_right.begin(), blocked_right.end(), 0);
      
      // search all unmatched nodes
      for (int i=0; i<n; i++) {
        if (edges_left[i]<0) {
          
          // attempt to find augmented path
          vector<int> path(1,-1);
          augment_left(path, i, m, edges_right, blocked_left, blocked_right);
          
          // if successful then augment
          if (path[0]>=0) {
            continue_augmenting = true;
            for (int j=0; j<int(path.size()/2); j++) {
              edges_left[path[j*2+1]] = path[j*2];
              edges_right[path[j*2]] = path[j*2+1];
            }
            number_assigned ++;
            
            // if best matching found then finish
            if (number_assigned==n) {
              return edges_left;
            }
          }
        }
      } // close i loop over unmatched nodes
    } // close while(continue_augmenting) loop
    
    // find minimum value in cost matrix, looking at all elements in which
    // neither the row nor the column is part of the minimum vertex cover
    double min_val = OVERFLO;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (blocked_left[i]==1 && blocked_right[j]==0 && m[i][j]<min_val) {
          min_val = m[i][j];
        }
      }
    }
    
    // add or subtract this value from cost matrix as required
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (blocked_left[i]==1 && blocked_right[j]==0) {
          m[i][j] -= min_val;
        }
        if (blocked_left[i]==0 && blocked_right[j]==1) {
          m[i][j] += min_val;
        }
      }
    }
    // at this point we have a new cost matrix and can repeat the process from the top
    
  } // rep loop
  
  // if reached this point then not managed to find best matching within given max_reps
  Rcpp::stop("Hungarian algorithm unable to find best matching"); // # nocov
  return edges_left;  // # nocov
}
