#ifndef __LINEAR_PROBLEM_SOLVER__
#define __LINEAR_PROBLEM_SOLVER__

#include "globaldefinitions.h"

bool solveBustingBoundaryPoints(dMatrix const &A, dVector const &b, dVector const &c, dVector &optSolution);
bool solveSimplexSyntheticBasis(dMatrix const &A, dVector const &b, dVector const &c, dVector &optSolution);
bool solveSimplex(dMatrix const &A, dVector const &b, dVector const &c, const dVector &x0, std::vector<int> &basisIndexes, dVector &optSolution);
void BuildDual(dMatrix const &A, dVector const &b, dVector const &c, dMatrix &dualMatrix, dVector &dualRightPart, dVector &dualC);

// everything less than this we will assume to be 0. 
const double machine_zero = 1e-13; 

#endif