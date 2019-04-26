%module pyUniDOE_swig
%{
#define SWIG_FILE_WITH_INIT
#include "wrapper.h"
%}

%include <std_vector.i>
%include <std_string.i>
%include <typemaps.i>

/* SWIG template for get_rand_list(int length) C++ routine */
namespace std {
  %template(DoubleVector) vector<double>;
  %template(DoubleVecVec) vector<vector<double> >;
  %template(IntegerVector) vector<int>;
  %template(IntegerVecVec) vector<vector<int> >;

}

%include "wrapper.h"
