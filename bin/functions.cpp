
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



// vector<double> compute_V(vd_t Int)
// {
//   vector<double> V(5);
//
//   double h1[5] = {  4.0, -4.0, -2.0,  2.0, 0.0};
//   double h2[5] = { 16.0, 16.0,  4.0,  4.0, 0.0};
//   double h3[5] = {  1.0, -1.0,  0.5, -0.5, 0.0};
//   double h4[5];
//   for (int i=0; i<5, i++){
//     h4[i] = (h2[i]/4.0-1.0)/3.0;
//
//
//     V[i] = VV00 + VV10 + csw*VV11 + csw*csw*VV12
//   }
//
//
//   return V;
// }
// double compute_S(vd_t Int)
// {
//
// }
