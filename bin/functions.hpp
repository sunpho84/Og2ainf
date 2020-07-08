#ifndef FUNC_HPP
#define FUNC_HPP


using namespace Eigen;
using namespace std;

double prod(Vector4d a);
double prod(Vector4d a,Vector4d b);
double prod(Vector4d a,Vector4d b,Vector4d c);
double prod(Vector4d a,Vector4d b,Vector4d c,Vector4d d);
double prod(Vector4d a,Vector4d b,Vector4d c,Vector4d d,Vector4d e);

Vector4d cwiseProduct(Vector4d a, Vector4d b);
Vector4d cwiseProduct(Vector4d a, Vector4d b, Vector4d c);
Vector4d cwiseProduct(Vector4d a, Vector4d b, Vector4d c, Vector4d d);
Vector4d cwiseProduct(Vector4d a, Vector4d b, Vector4d c, Vector4d d, Vector4d e);


Matrix4d make_propagator(Vector4d kt2_dir, Vector4d kt4_dir, Vector4d kt6_dir);

#endif
