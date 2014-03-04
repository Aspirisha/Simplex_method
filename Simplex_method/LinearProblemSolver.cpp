#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "linearproblemsolver.h"
#include "ioInterface.h"
#include "gaussolver.h"

using namespace boost::numeric::ublas;

static void _getBasisCombinations(int start, int m, int n, std::vector<int> &chosen, std::vector<std::vector<int> > &variants);

std::vector<std::vector<int> > getBasisCombinations(int m, int n)
{
  std::vector<std::vector<int> > retVal;
  std::vector<int> v(n);

  _getBasisCombinations(0, m, n, v, retVal);
  return retVal;
}

// solves canonical linear programming problem bysting all the boundary points
// problem: 
// Ax = b
// c * x -> min
bool solveBustingBoundaryPoints(dMatrix const &A, dVector const &b, dVector const &c, dVector &optSolution)
{
  int m = A.size1();
  int n = A.size2();

  std::vector<std::vector<int> > basisCombinations = getBasisCombinations(m, n);
  int numberOfcombinations = basisCombinations.size();
  dMatrix cur(m, m);
  double optimalValue = 1e307; // very big to init
  
  bool opt_found = false;
  for (int i = 0; i < numberOfcombinations; i++)
  {
    int k = 0;
    for (int j = 0; j < m; j++)
    {
      while (basisCombinations[i][k] == 0)
        k++;
      column(cur, j) = column(A, k); 
      k++;
    }
    dVector x_trancated = gausSolve(cur, b);
    if (x_trancated.size() != m) // matrix cur is with 0 det
      continue;

    dVector x = zero_vector<double>(n);
    k = 0;

    bool negative_component_detected = false;
    for (int j = 0; j < n; j++)
    {
      if (basisCombinations[i][j] == 1)
      {
        x(j) = x_trancated(k++);
        if (x(j) < 0)
          negative_component_detected = true;
      }
    }
    if (negative_component_detected)
      continue;
    double currentValue = inner_prod(x, c);
    if (currentValue < optimalValue)
    {
      optimalValue = currentValue;
      optSolution = x;
      opt_found = true;
    }
  }
  return opt_found;
}

bool solveSimplexSyntheticBasis(dMatrix const &A, dVector const &b, dVector const &c, dVector &optSolution)
{
  int m = A.size1();
  int n = A.size2();
  int new_n = n + m;

  dVector x0 = zero_vector<double>(new_n);
  dVector new_c = zero_vector<double>(new_n);
  dMatrix B = zero_matrix<double>(m, new_n);

  for (int i = 0; i < n; i++)
    column(B, i) = column(A, i);

  std::vector<int> basisIndexes(m);
  for (int i = n; i < new_n; i++)
  {
    B(i - n, i) = 1;
    x0(i) = b(i - n);
    basisIndexes[i - n] = i;
    new_c(i) = 1; // c == (0, 0, ..., 0, 1, 1, ... 1)
  }
  
  dVector startPoint(new_n);
  __DUMP__VARIABLE__(x0, "startPoint");
  solveSimplex(B, b, new_c, x0, basisIndexes, startPoint);
  double mu = inner_prod(new_c, startPoint);
  __DUMP__VARIABLE__(startPoint, "startPoint");
  __DUMP__VARIABLE__(startPoint, "startPoint");
  if (mu > 1e-13) // no solution
    return false;
  startPoint.resize(n);
  __DUMP__VARIABLE__(startPoint, "x0");
  solveSimplex(A, b, c, startPoint, basisIndexes, optSolution);
  return true;
}

// x0 is start boundary point
bool solveSimplex(dMatrix const &A, dVector const &b, dVector const &c, const dVector &x0, std::vector<int> &basisIndexes, dVector &optSolution)
{
  bool solutionFound = false;
  bool solution_exists = false;
  int m = A.size1();
  int n = A.size2();

  __DUMP__VARIABLE__(A, "A");
  __DUMP__VARIABLE__(b, "b");
  dVector x = x0; // current bounder point. 
  dVector cb(m);
  dMatrix B = zero_matrix<double>(m, m); // square matrix built from basis columns of A for current boundary point
  std::vector<double> delta(n);

  while (!solutionFound)
  {
    int j = 0;
    
    for (int i = 0; i < m; i++)
    {
      int ind = basisIndexes[i];
      column(B, j) = column(A, ind); // making square m x m matrix using indexes from basisIndex array
      cb(j) = c(ind); // taking appropriate c(j)
      j++;
    }
    
    __DUMP__VARIABLE__(B, "B");
    dMatrix invB = gausInv(B);
    __DUMP__VARIABLE__(invB, "invB");

    bool allNotPositive = true;
    int positive_k_index = -1;
    dVector positive_x1;
    for (int i = 0; i < n; i++)
    {
      dVector x1 = prod(invB, column(A, i));
      __DUMP__VARIABLE__(x1, "x1");
      delta[i] = inner_prod(cb, x1) - c(i);
      // TODO manage with machine 0:
      if (delta[i] > 1e-10) // it automatically means that i is not basis index, cause else delta[i] = 0. change 1e-10 to machine_epse then
      {
        allNotPositive = false;
        bool x1NotPositive = true;
        for (int k = 0; k < m; k++) // test if there is component in invB * A[j]
          if (x1(k) > 0)
            x1NotPositive = false;
        if (x1NotPositive) // <c, x> is not bounded cause invB * A[j] <= 0 for j - non basis index
          return 1;

        if (positive_k_index == -1)
        {
          positive_k_index = i;
          positive_x1 = x1;
        }
      }
    }
    __DUMP__VARIABLE__(positive_x1, "positive_x1");
    __DUMP__VARIABLE__(positive_k_index, "positive component: ");
    //If all deltas are <= 0 then x - opt, exit from loop
    if (allNotPositive)
      solutionFound = true;
    
    // need to make next iteration, no solution yet was found
    if (!solutionFound)
    {
      int excludingIndex = -1;
      std::vector<int> I; // building Ik - indexes of x1_positive which coords are > 0
      for (int i = 0; i < m; i++)
      {
        if (positive_x1(i) > 0)
          I.push_back(i);
      }

      int iter = -1;
      dVector x3;
      double lambda;
      double new_component = 0;
      while (I.size() >= 1)
      {
        std::vector<int> new_I;
        lambda = 1e307;

        if (iter == -1)
          x3 = prod(invB, b);
        else
          x3 = prod(invB, column(B, iter));
        iter++;
        __DUMP__VARIABLE__(x3, "x3");
        
        for (size_t i = 0; i < I.size(); i++)
        {
          int j = I[i];
          double cur_lambda = x3(j) / positive_x1(j);
          if (lambda > cur_lambda)
          {
            new_I.clear();
            lambda = cur_lambda;
            new_I.push_back(j);
          } 
          else if (lambda == cur_lambda)
            new_I.push_back(j);
        }
        
        if (iter == 0)
          new_component = lambda;
        I = new_I;
        if (I.size() == 1) // point is not degenerate so stop this process
          break;
      }

      int index_to_delete = I[0];  
      
      //build new bounder point xk
      dVector xk = zero_vector<double>(n);
      __DUMP__VARIABLE__(x, "x");
     
      for (int i = 0; i < m; i++)
      {
        int j = basisIndexes[i];
        xk(j) = x(j) - new_component * positive_x1(i);
      }
      // insert new basis index instead of old one (I[0] is old, positive_k_index is new)
      basisIndexes[index_to_delete] = positive_k_index;

      xk(positive_k_index) = new_component;
      __DUMP__VARIABLE__(xk, "xk");
      x = xk;
      solution_exists = true;
    }
  }
  optSolution = x;
  return solution_exists; 
}



void _getBasisCombinations(int start, int m, int n, std::vector<int> &chosen, std::vector<std::vector<int> > &variants)
{
  if (m == 0)
  {
    variants.push_back(chosen);
    return;
  }

  for (int i = start; i < n; i++)
  {
    if (!chosen[i])
    {
      chosen[i] = 1;
      _getBasisCombinations(i + 1, m - 1, n, chosen, variants);
      chosen[i] = 0;
    }
  }
}

/**
 * @param A is matrix of initial linear programming problem. It's assumed that problem is represented canonically.
 * @param b is vector representing right part of problem.
 * @param c is parameter of minimization
 * @param dualMatrix is output matrix (matrix for dual problem)
 * @apram dualRightPart --||--
 * @param dualC --||--
 */ 
void BuildDual(dMatrix const &A, dVector const &b, dVector const &c, dMatrix &dualMatrix, dVector &dualRightPart, dVector &dualC)
{
  dualMatrix = -trans(A);
  dualRightPart = -c;
  dualC = -b;
}