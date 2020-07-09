
#include "global.hpp"

using namespace std;
using namespace Eigen;

/* Volume */
int L;
int T;
Vector2d dim;
Vector4d V;
int APBC;
int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

/* Action */
double c1;
double r;
double csw;
double alpha;
double beta;
int action;

/* External momentum */
Vector4d np;
Vector4d ap;
double Np0;

/* Loop momentum */
double al;

/* Invariants */
Vector4d Id(1.0,1.0,1.0,1.0);
