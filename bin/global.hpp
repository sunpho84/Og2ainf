#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include <Eigen/Dense>
#include <valarray>
#include <complex>
#include <array>
#include <array>
#include <vector>

using namespace std;
using namespace Eigen;

/* Volume */
extern int L;
extern int T;
extern array<int,2> dim;
extern array<int,4> V;
extern int APBC;
extern int kron_delta[4][4];

/* Action */
extern double c1;
extern double r;
extern double csw;
extern double alpha;
extern double beta;
extern int action;

/* External momentum */
extern array<double,4> np;
extern array<double,4> ap;
extern double Np0;

/* Loop momentum */
extern double al;

/* Invariants */
// extern Vector4d Id;


#endif
