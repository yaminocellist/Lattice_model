#pragma once

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>
#include <valarray>
#include <vector>
#include <complex>
#include <algorithm>

#include "settings.h"
// #include </opt/homebrew/Cellar/boost/1.78.0_1/include/boost/math/special_functions/bessel.hpp>
// #include </opt/homebrew/Cellar/boost/1.78.0_1/include/boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>




void printvd (vector<double> vd);
void printvd (vector<int> vd);
void print2dv (vector2d v);
void print3dv (vector3d v);
int first_non_zero (vector<double> v);
vector<MatrixXd> set_Qn_delta (Parameters p, int n, int g);
#endif
