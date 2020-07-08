
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "global.hpp"

using namespace Eigen;
using namespace std;


// double k(const int& i, const int& j, const vector<int>& dim)
// {
//   return 2.0*M_PI/dim[i]*(j+0.5) - al*sin(2.0*M_PI/dim[i]*(j+0.5));
// }
// double kn(const int i, const int j, const vector<int>& dim, const vector<double>& ap, )
// {
//   return k(i,j,dim) + 2.0*ap[i];
// }
// double kp(const int i, const int j)
// {
//   return k(i,j) + ap[i];
// }
// double J(const int& i, const int& j)
// {
//   return 1.0 - al*cos(2.0*M_PI/dim[i]*(j+0.5));
// }
