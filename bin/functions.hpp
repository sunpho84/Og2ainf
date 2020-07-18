#ifndef FUNC_HPP
#define FUNC_HPP

#include "global.hpp"
#include "simd.hpp"
#include <array>

// using namespace Eigen;
using namespace std;

template <typename T>
inline array<array<T,4>,4> make_propagator(array<T,4> kt2_dir, array<T,4> kt4_dir, array<T,4> kt6_dir)
{
  array<array<T,4>,4> ktpo2, ktso2, A;

  T kt2=T(), kt4=T(), kt6=T();

  for(int mu=0; mu<4; mu++)
  {
    kt2 += kt2_dir[mu];
    kt4 += kt4_dir[mu];
    kt6 += kt6_dir[mu];
  }

  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++){

      for(int x=0;x<N;x++)
      {
        ktpo2[mu][nu][x]=1.0;
        ktso2[mu][nu][x]=0.0;
      }

      for(int rho=0;rho<4;rho++){
        if(mu!=rho && nu!=rho)
        {
          ktpo2[mu][nu]*=kt2_dir[rho];
          ktso2[mu][nu]+=kt2_dir[rho];
        }
      }
    }

  T kt22=kt2*kt2;
  T kt23=kt2*kt2*kt2;
  T kt42=kt4*kt4;

  T Deltakt=(kt2-c1*kt4)*(kt2-c1*(kt22+kt4)+(Real)0.5*c12*(kt23+(Real)2.0*kt6-kt2*kt4));
  for(int rho=0;rho<4;rho++) Deltakt-=(Real)4.0*c13*kt4_dir[rho]*ktpo2[rho][rho];

  for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++){
	    A[mu][nu]=((Real)1.0 - kron_delta[mu][nu])/Deltakt*(kt22-c1*kt2*((Real)2.0*kt4+kt2*ktso2[mu][nu])+c12*(kt42+kt2*kt4*ktso2[mu][nu]+kt22*ktpo2[mu][nu]));
	    A[mu][nu]-=((Real)1.0 - kron_delta[mu][nu]);
    }

  return A;
}

vector<Real> compute_Gamma(Real *Int, Real *IntS);
vector<Real> compute_Z(Real *Int, Real *IntS);

inline Real norm4(const array<Real,4>& mom)
{
  return mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]+mom[3]*mom[3];
}

void find_eqmoms();

#endif
