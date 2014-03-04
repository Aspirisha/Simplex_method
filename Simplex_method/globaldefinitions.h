#ifndef __GLOBAL_DEFINITIONS__H__
#define __GLOBAL_DEFINITIONS__H__

// uncomment line below to see debug dumps
#define __DEBUG__MODE__ON__ 

#include <vector>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

typedef boost::numeric::ublas::matrix<double> dMatrix; // matrix of doubles
typedef boost::numeric::ublas::vector<double> dVector;
typedef std::vector<double> line;

#endif