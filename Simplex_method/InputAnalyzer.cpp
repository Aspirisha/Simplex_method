#include "inputanalyzer.h"
#include <iostream>
#include <map>
#include <string>
#include <boost/regex.hpp>

using namespace std;
using namespace boost::numeric::ublas;

typedef enum { LESS_EQ, GREATER_EQ, EQ } sign_types;

static void makeCanonical(dMatrix &A, dVector &b, std::vector<string> &varNames, int &extra_variables_number, const std::vector<sign_types> &sign_type);

bool parseLE(const std::vector<string> &problemText, dMatrix &A, dVector &b, std::vector<string> &varNames, int &extra_variables_number)
{
  map<string, int> existingVariables; // map containing variables which were already seen in problem text by now. For each var it contains it's number in order they appear
  int variablesAmount = 0;

  int m = problemText.size();
  std::vector<std::vector<double> > x(1);
  std::vector<sign_types> sign_type; // 1: <=, 2: >=, 3: =
  std::vector<double> c(1);
  int trueLinesNumber = 0;
  bool isInvisiblePlus = true; // first plus in an expression can be skipped, it's not an error
  int pos = 0;  // position of regexp pointer in current string
  double sign = 1; // sign will change to -1 after we pass the equals-symbol
  bool equalsFound = false;
  char operation = '+';

  string regexpStr = "\\s*([+-]?)\\s*(\\d+[.,]\\d+|\\d*)\\s*([a-z]*[0-9]*)\\s*(=|>=|<=)?([;\\n]?)";
  boost::smatch what;
  boost::regex regexp(regexpStr);
 
  for (int lineNumber = 0; lineNumber < m; lineNumber++)
  {
    std::string::const_iterator start = problemText[lineNumber].begin();
    std::string::const_iterator end   = problemText[lineNumber].end();

    bool endOfLineFound = false;
    while (boost::regex_search(start, end, what, regexp))
    {
      if (what[3] == "" && what[2] == "") // no matches indeed
        break;

      start = what[5].second;

      int varNumber = -1; // number of captured variable in existingVariables array. 
      string varName(what[3]);

      if (varName != "") 
      {
        if (!existingVariables.count(varName)) 
        {
            existingVariables[varName] = variablesAmount;
            for (int i = 0; i <= lineNumber; i++)
              x[i].push_back(0);
            varNumber = variablesAmount;
            ++variablesAmount;
        } 
        else 
          varNumber = existingVariables[varName];
      }

      string current_sign(what[1]);
      string current_coef(what[2]);
      if (current_sign == "+" || current_sign == "-")
        operation = current_sign[0];
      else if (isInvisiblePlus)
      {
        operation = '+';
      }
      else
        return false;

      isInvisiblePlus = false; // no more invisible pluses unless we find right part of equation

      double cur_summand = sign;
      if (operation == '-')
        cur_summand *= -1;

      if (current_coef != "")
        cur_summand *= atof(current_coef.c_str());
      else if (varName == "") // no coefficient and no varname - it's error
        return false;

      if (varName == "") 
        c[lineNumber] += cur_summand;
      else 
        x[lineNumber][varNumber] += cur_summand;
         

      string sign_of_expression(what[4]);
      if (sign_of_expression != "") 
      { // '=' or '>=' or '<='
        if (sign_of_expression == "=")
          sign_type.push_back(EQ);
        else if (sign_of_expression == ">=")
          sign_type.push_back(GREATER_EQ);
        else
          sign_type.push_back(LESS_EQ);

        sign = -1;
        isInvisiblePlus = true;
        if (equalsFound)
          return false; // more then one equals sign is found
        equalsFound = true;
      }
        
      string endofLine(what[5]);
      if (endofLine != "") // ';' or '\n'
      {
        if (isInvisiblePlus || !equalsFound) // impossible unless right or left part of equation is empty, then it's error
          return false;
        x.push_back(std::vector<double>(variablesAmount));
        c.push_back(0);
        sign = 1;
        isInvisiblePlus = true;
        equalsFound = false;
        endOfLineFound = true;
      }
    }

    if (!endOfLineFound)
      return false;
  }

  A.resize(m, variablesAmount);
  b.resize(m);

  for (int i = 0; i < m; i++) 
    b(i) = -c[i];
  
  for (int i = 0; i < m; i++) 
    for (int j = 0; j < variablesAmount; j++) 
      A(i, j) = x[i][j];

 
  // facepalm, but it was fast to code XD. TODO avoid O(n^2)
  for (int i = 0; i < existingVariables.size(); i++) 
     for (map<string,int>::iterator it = existingVariables.begin(); it != existingVariables.end(); ++it) 
       if (it->second == i)
         varNames.push_back(it->first);

  // now it's time to make problem canonical
  makeCanonical(A, b, varNames, extra_variables_number, sign_type);
  
  return true;
}

void makeCanonical(dMatrix &A, dVector &b, std::vector<string> &varNames, int &extra_variables_number, const std::vector<sign_types> &sign_type)
{
  extra_variables_number = 0;
  int m = A.size1();
  int n = A.size2();

  int eq_number = m;
  for (int i = 0; i < eq_number; i++)
  {
    if (sign_type[i] == LESS_EQ) // change sign of line to >= (simply multiply by -1 all coefficients)
    {
      for (int j = 0; j < n; j++)
        A(i, j) = -A(i, j);
      b(i) = -b(i);
    }

    if (b(i) < 0) // add extra variable vi = -b[i]
    {
      A.resize(++m, ++n);
      for (int j = 0; j < n - 1; j++)
        A(m - 1, j) = 0;
      A(m - 1, n - 1) = 1;
      for (int j = 0; j < m - 1; j++)
        A(j, n - 1) = 0;
      A(i, n - 1) = 1;
      b.resize(m);
      b(m - 1) = -b(i);
      b(i) = 0;
      b(i) = 0;

      std::stringstream out;
      out << extra_variables_number++ << "u"; // starting extra variable name from number guarantees it was not reserved by user
      string varName = out.str();
      varNames.push_back(varName);
    }

    if (sign_type[i] == EQ)
      continue;
    
    //create new variable and alter matrix
    std::stringstream out;
    out << extra_variables_number++ << "u"; // starting extra variable name from number guarantees it was not reserved by user
    string varName = out.str();
    varNames.push_back(varName);

    A.resize(m, ++n);
    for (int j = 0; j < m; j++)
      A(j, n - 1) = 0;
    A(i, n - 1) = -1;
  }
}