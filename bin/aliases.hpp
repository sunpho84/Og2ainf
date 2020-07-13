#ifndef ALIASES_HPP
#define ALIASES_HPP

#ifdef HAVE_CONFIG_H
#include <config.hpp>
#endif

// #ifndef EIGEN_H
// #include <Eigen/Dense>
// #define EIGEN_H
// #endif

#include <valarray>
#include <complex>
#include <array>
#include <vector>
#include <tuple>

using namespace std;
//using namespace Eigen;

// //coordinates in the lattice
// using coords_t=array<int,4>;
//
// //p components in the lattice
// using p_t=array<double,4>;
//
// //complex double
// using dcompl=complex<double>;
//
// //propagator (12X12)
// using prop_t=Matrix<dcompl,12,12>;
//
// //list of propagators
// using vprop_t=valarray<prop_t>;
// using vvprop_t=valarray< vprop_t >;
// using vvvprop_t=valarray< vvprop_t >;
// using vvvvprop_t = valarray< vvvprop_t >;
//
// //list of gamma for a given momentum
// using qline_t=valarray<prop_t>;
// using vqline_t=valarray<qline_t>;
// using vvqline_t=valarray<vqline_t>;
// using vert_t = vvqline_t;
//
// //list of jackknife propagators
// using jprop_t=vvvprop_t;
//
// //list of jackknife vertices
// using jvert_t=valarray< vert_t > ;

// //valarray of complex double
// using vd_t=valarray<double>;
// using v1d_t=vd_t;
//
// //valarray of valarray of complex double
// using vvd_t=valarray< valarray<double> > ;
// using v2d_t=vvd_t;
//
// //valarray of valarray of valarray of complex double
// using vvvd_t=valarray< valarray< valarray<double> > >;
// using vvvvd_t=valarray<vvvd_t>;
// using vvvvvd_t=valarray<vvvvd_t>;
// using vvvvvvd_t=valarray<vvvvvd_t>;
// using v3d_t=vvvd_t;
// using v4d_t=vvvvd_t;
// using v5d_t=vvvvvd_t;
// using v6d_t=vvvvvvd_t;


#endif
