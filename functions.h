#pragma once

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>
#include <valarray>
#include <vector>
#include <complex>
#include <algorithm>

#include "settings.h"
#include </opt/homebrew/Cellar/boost/1.78.0_1/include/boost/math/special_functions/bessel.hpp>
#include </opt/homebrew/Cellar/boost/1.78.0_1/include/boost/math/quadrature/gauss_kronrod.hpp>

using std::vector;
using vector2d = vector<vector <double>>;

void printvd (vector<double> vd);
void printvd (vector<int> vd);
void print2dv (vector2d v);
int sum_mg (vector<int> mg, int it);
int find_g_1 (int i, vector<int> mg, int f);
int first_non_zero (vector<double> v);
vector<int> rest_mg (vector<int> mg, int g);
#endif