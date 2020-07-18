#ifndef GLOBAL_HPP
#define GLOBAL_HPP

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
extern double c1;
extern double c12;
extern double c13;
extern double r;
extern double csw;
extern double alpha;
extern double beta;
extern int action;

/* External momentum */
// extern array<double,4> np;
// extern array<double,4> ap;
extern double Np0;

extern vector<array<int,4>> np_list;
extern vector<array<double,4>> ap_list;
extern vector<array<double,4>> ap_eq;
extern vector<int> tag_list;
extern int moms;
extern int eqmoms;

/* Loop momentum */
extern double al;

/* Invariants */
// extern Vector4d Id;


#endif
