#ifndef _SIMPLEX_INPUT_ANALYZER_H__
#define _SIMPLEX_INPUT_ANALYZER_H__

#include "globaldefinitions.h"

bool parseLE(const std::vector<std::string> &problemText, dMatrix &A, dVector &b, std::vector<std::string> &varNames, int &extra_variables_number);

#endif