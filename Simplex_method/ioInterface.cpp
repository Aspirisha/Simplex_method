#include "ioInterface.h"
#include "inputanalyzer.h"
#include <fstream>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

void __DUMP__VARIABLE__(const dMatrix &A, const string &name)
{
  #ifdef __DEBUG__MODE__ON__
  cout << name << " = " << A << endl;
  #endif
}

void __DUMP__VARIABLE__(const dVector &b, const string &name)
{
  #ifdef __DEBUG__MODE__ON__
  cout << name << " = " << b << endl;
  #endif
}

void __DUMP__VARIABLE__(int n, const string &name)
{
  #ifdef __DEBUG__MODE__ON__
  cout << name << " = " << n << endl;
  #endif
}

bool readConditions(dVector &c, dMatrix &A, dVector &b, int &extra_variables_number, std::vector<string> &varNames)
{
  int max_dimension = 10;
  int max_condition_number = 10;
  
  cout << "Insert vector c:\n";
  cout << ">> ";

  string s;
  getline(cin, s);
  stringstream ss(s);
  int n = 0;
  c.resize(max_dimension);
  double d;
  while ((ss >> d) && (n < max_dimension))
    c(n++) =  d;
  
  if (n < max_dimension)
    c.resize(n);

  int m = 0; // restrictions number
  s = "Not empty";
  std::vector<string> conditions;

  cout << "Insert linear constraints:\n";
  cin.sync();
  while (m < max_condition_number)
  {
    cout << ">> ";
    getline(cin, s);
    if (s.empty())
      break;
    conditions.push_back(s);
    m++;
  }
  
  if (!parseLE(conditions, A, b, varNames, extra_variables_number))
    return false;
  c.resize(n + extra_variables_number);
  if (c.size() != A.size2())
    return false;
  for (int i = 0; i < extra_variables_number; i++)
    c(n + i) = 0;
  __DUMP__VARIABLE__( A, "A");
  __DUMP__VARIABLE__( b, "b");
  return true;
}

void printHelp()
{
  ifstream in("help.txt");
  string s;
  while (getline(in, s))
    cout << s << endl;
  cout << endl;
  in.close();
}

void printAnswer(const dVector &x, const std::vector<string> &varNames, const string &methodName)
{
  cout << "Solution via " << methodName <<":\n";
  int n = x.size();
  for (int i = 0; i < n; i++)
    cout << varNames[i] << " = " << x(i) << endl;
}

void printError(const std::string &methodName)
{
  cout << "No solution was found or error has occured while solving problem with " << methodName << "\n";
}