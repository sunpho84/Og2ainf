#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include "simd.hpp"
// #include <Eigen/Dense>
#include <valarray>
#include <complex>
#include <array>
#include <vector>

using namespace std;
// using namespace Eigen;

/* Volume */
extern int L;
extern int T;
extern array<int,2> dim;
extern array<int,4> V;
extern int APBC;
extern int kron_delta[4][4];

extern int eq;

/* Action */
extern Real c1;
extern Real c12;
extern Real c13;
extern Real r;
extern Real csw;
extern Real alpha;
extern Real beta;
extern int action;

/* External momentum */
// extern array<Real,4> np;
// extern array<Real,4> ap;
extern Real Np0;

extern vector<array<int,4>> np_list;
extern vector<array<Real,4>> ap_list;
extern vector<array<Real,4>> ap_eq;
extern vector<int> tag_list;
extern int moms;
extern int eqmoms;

/* Loop momentum */
extern Real al;

/* Invariants */
// extern Vector4d Id;


#endif
