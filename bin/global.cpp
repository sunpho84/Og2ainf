#include "global.hpp"

using namespace std;
// using namespace Eigen;

/* Volume */
int L;
int T;
array<int,2> dim;
array<int,4> V;
int APBC;
int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

int eq;

/* Action */
Real c1;
Real c12;
Real c13;
Real r;
Real csw;
Real alpha;
Real beta;
int action;

/* External momentum */
// array<Real,4> np;
// array<Real,4> ap;
Real Np0;

vector<array<int,4>> np_list;
vector<array<Real,4>> ap_list;
vector<array<Real,4>> ap_eq;
vector<int> tag_list;
int moms;
int eqmoms;

/* Loop momentum */
Real al;

/* Invariants */
// Vector4d Id(1.0,1.0,1.0,1.0);
