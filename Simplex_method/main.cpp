#include <iostream>

#include "iointerface.h"
#include "linearproblemsolver.h"

using namespace std;

int main(void)
{
  dMatrix A;
  dVector b;
  dVector c;
  int extra_variables_number = 0;
  std::vector<string> varNames;

  printHelp();
  if (!readConditions(c, A, b, extra_variables_number, varNames))
  {
    cout << "Wrong input\n";
    return 0;
  }

  int m = A.size1(); // number of rows
  int n = A.size2(); // number of columns

  dVector bustingSolution(n);
  if (!solveBustingBoundaryPoints(A, b, c, bustingSolution))
    printError("busting boundary points method");
  else
  {
    bustingSolution.resize(n - extra_variables_number);
    printAnswer(bustingSolution, varNames, "busting boundary points method");
  }

  dVector simplexSolution(n);
  if (!solveSimplexSyntheticBasis(A, b, c, simplexSolution))
    printError("simplex method");
  else
  {
    simplexSolution.resize(n - extra_variables_number);
    printAnswer(bustingSolution, varNames, "simplex method");
  }
  return 0;
}