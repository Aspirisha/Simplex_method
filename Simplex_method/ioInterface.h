#ifndef __IO_INTERFACE_H__
#define __IO_INTERFACE_H__

#include "globaldefinitions.h"

void __DUMP__VARIABLE__(const dMatrix &A, const std::string &name);
void __DUMP__VARIABLE__(const dVector &b, const std::string &name);
void __DUMP__VARIABLE__(int n, const std::string &name);
bool readConditions(dVector &c, dMatrix &A, dVector &b, int &extra_variables_number, std::vector<std::string> &varNames);
void printHelp();
void printAnswer(const dVector &x, const std::vector<std::string> &varNames, const std::string &methodName);
void printError(const std::string &methodName);
#endif