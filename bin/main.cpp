#ifdef HAVE_CONFIG_H
#include <config.hpp>
#endif

#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <array>
#include <vector>
#include <numeric>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
// #include <Eigen/Dense>
// #include <Eigen/StdVector>
#include "aliases.hpp"
#include "global.hpp"
#include "read_input.hpp"
#include "functions.hpp"

#include "simd.hpp"

#define EulerGamma 0.5772156649015329
#define F0 4.36923
#define Z0 0.154933

array<Real,4> p;
Real P2;
Real PxSINP;
vector<Real> j,k,kt,kt2_v,kt4_v,kt6_v,s_v,ss_v,s2_v,ss2_v,c_v,cc_v,
kn,kp,sp_v,cp_v,sp2_v,sq2_v,ssp_v,ssp2_v;

template <typename F>
void intComp(std::array<F,12>& Int,std::array<F,10>& IntS,std::array<int,N> i)
{
  std::array<int,N> i4,i3,i2,i1;
  for(int j=0;j<N;j++)
  {
    i4[j]=i[j]%T;
    i[j]/=T;
    i3[j]=i[j]%L;
    i[j]/=L;
    i2[j]=i[j]%L;
    i[j]/=L;
    i1[j]=i[j];
  }

  array<F,4> kt2, kt4, kt6;
  array<F,4> s,ss,s2,ss2,c,cc;
  array<F,4> sp,cp,sp2,sq2,ssp,ssp2;
  F S2=F(), SP2=F(), SQ2=F(), SS2=F(), SSP2=F();
  F Fac;
  F Delta1, Delta2, Delta1p, Delta2p;
  F Delta12, Delta13, Delta14, Delta22, Delta1p2, Delta2p2, Delta12p2, Den;
  F SxSP=F(), CxCP=F(), SSxSSP=F(), S2xSP2=F(), S2xSS2=F(), S2xSSP2=F(), SP2xSSP2=F(), S2xSSxSSP=F(), CPxSxSS=F(), CPxSxSSP=F(), CPxSPxSSP=F(), CxCPxSSxSSP=F(), CxSPxSSP=F(), CCxS2xSS2=F();
  F PxSS=F(), PxSSP=F(), PxCPxS=F(), PxS2xSSP=F(), PxS2xSS=F(), PxCxCPxSS=F(), PxCPxSP=F(), PxCxSP=F(), PxSP2xSSP=F();
  F CCxS2xSS2xAtxS2=F(), S2xSS2xAtxSS2=F(), S2xSS2xAtxS2=F(), SS2xAtxSS2=F(), SS2xAtxS2=F(), S2xAtxS2=F(), AtxS2=F();
  F CPxSxSSPxAtxCPxSxSSP=F(), CPxSxSSxAtxSSxSSP=F(), CPxSxSSPxAtxSxSP=F(), CPxSPxSSPxAtxS2=F(), CxCPxSSxSSPxAtxS2=F(), CxSPxSSPxAtxS2=F(), CxCPxAtxS2=F(), S2xSSxSSPxAtxS2=F(), S2xSSP2xAtxS2=F(), SP2xSSP2xAtxS2=F(), SSxSSPxAtxSSxSSP=F(), SSxSSPxAtxSxSP=F(), SSxSSPxAtxSS2=F(), SSxSSPxAtxS2=F(), SxSPxAtxSxSP=F(), SSP2xAtxS2=F(), SP2xAtxS2=F();
  F PxCPxSxAtxCPxSxSSP=F(), PxCPxSxAtxSxSP=F(), PxCPxSPxAtxS2=F(), PxCxCPxSSxAtxS2=F(), PxCxSPxAtxS2=F(), PxSSPxAtxS2=F(), PxSP2xSSPxAtxS2=F(), PxS2xSSPxAtxS2=F(), PxS2xSSxAtxS2=F(), PxSSxAtxCPxSxSS=F(), PxSSxAtxSSxSSP=F(), PxSSxAtxSxSP=F(), PxSSxAtxSS2=F(), PxSSxAtxS2=F();
  array<array<F,4>,4> A;

  for(int x=0;x<N;x++)
  {
    Fac[x] = j[i1[x]]*j[i2[x]]*j[i3[x]]*j[i4[x]+L]*16.0*M_PI*M_PI/((Real)(L*L*L*T));

    auto fill=[x,i1,i2,i3,i4](array<F,4>& kt,const vector<Real>& kt_v)
    {
      kt[0][x] = kt_v[i1[x]];
      kt[1][x] = kt_v[i2[x]];
      kt[2][x] = kt_v[i3[x]];
      kt[3][x] = kt_v[i4[x]+L];
    };
    fill(kt2,kt2_v);
    fill(kt4,kt4_v);
    fill(kt6,kt6_v);
    fill(s,s_v);
    fill(ss,ss_v);
    fill(s2,s2_v);
    fill(ss2,ss2_v);
    fill(c,c_v);
    fill(cc,cc_v);

    auto fill3=[x,i1,i2,i3,i4](array<F,4>& t,const vector<Real>& t_v)
    {
      t[0][x] = t_v[i1[x]];
      t[1][x] = t_v[i2[x]+L];
      t[2][x] = t_v[i3[x]+2*L];
      t[3][x] = t_v[i4[x]+3*L];
    };
    fill3(sp,sp_v);
    fill3(cp,cp_v);
    fill3(sp2,sp2_v);
    fill3(sq2,sq2_v);
    fill3(ssp,ssp_v);
    fill3(ssp2,ssp2_v);

    A = make_propagator(kt2,kt4,kt6);
  }

  // asm("#CIAO!");

  // #pragma unroll
  for(int mu=0; mu<4; mu++){
    S2   += s2[mu];
    SP2  += sp2[mu];
    SQ2  += sq2[mu];
    SS2  += ss2[mu];
    SSP2 += ssp2[mu];

    SxSP        += s[mu]*sp[mu];
    CxCP        += c[mu]*cp[mu];
    SSxSSP      += ss[mu]*ssp[mu];
    S2xSP2      += s2[mu]*sp2[mu];
    S2xSS2      += s2[mu]*ss2[mu];
    S2xSSP2     += s2[mu]*ssp2[mu];
    SP2xSSP2    += sp2[mu]*ssp2[mu];
    S2xSSxSSP   += s2[mu]*ss[mu]*ssp[mu];
    CPxSxSS     += cp[mu]*s[mu]*ss[mu];
    CPxSxSSP    += cp[mu]*s[mu]*ssp[mu];
    CPxSPxSSP   += cp[mu]*sp[mu]*ssp[mu];
    CxCPxSSxSSP += c[mu]*cp[mu]*ss[mu]*ssp[mu];
    CxSPxSSP    += c[mu]*sp[mu]*ssp[mu];
    CCxS2xSS2   += cc[mu]*s2[mu]*ss2[mu];

    PxSS      += p[mu]*ss[mu];
    PxSSP     += p[mu]*ssp[mu];
    PxCPxS    += p[mu]*cp[mu]*s[mu];
    PxS2xSSP  += p[mu]*s2[mu]*ssp[mu];
    PxS2xSS   += p[mu]*s2[mu]*ss[mu];
    PxCxCPxSS += p[mu]*c[mu]*cp[mu]*ss[mu];
    PxCPxSP   += p[mu]*cp[mu]*sp[mu];
    PxCxSP    += p[mu]*c[mu]*sp[mu];
    PxSP2xSSP += p[mu]*sp2[mu]*ssp[mu];

    for(int nu=0; nu<4; nu++)
    {
      CCxS2xSS2xAtxS2 += cc[mu]*s2[mu]*ss2[mu]*A[mu][nu]*s2[nu];
      S2xSS2xAtxSS2   += s2[mu]*ss2[mu]*A[mu][nu]*ss2[nu];
      S2xSS2xAtxS2    += s2[mu]*ss2[mu]*A[mu][nu]*s2[nu];
      SS2xAtxSS2      += ss2[mu]*A[mu][nu]*ss2[nu];
      SS2xAtxS2       += ss2[mu]*A[mu][nu]*s2[nu];
      S2xAtxS2        += s2[mu]*A[mu][nu]*s2[nu];
      AtxS2           += A[mu][nu]*s2[nu];

      CPxSxSSPxAtxCPxSxSSP += cp[mu]*s[mu]*ssp[mu]*A[mu][nu]*cp[nu]*s[nu]*ssp[nu];
      CPxSxSSxAtxSSxSSP  += cp[mu]*s[mu]*ss[mu]*A[mu][nu]*ss[nu]*ssp[nu];
      CPxSxSSPxAtxSxSP += cp[mu]*s[mu]*ssp[mu]*A[mu][nu]*s[nu]*sp[nu];
      CPxSPxSSPxAtxS2  += cp[mu]*sp[mu]*ssp[mu]*A[mu][nu]*s2[nu];
      CxCPxSSxSSPxAtxS2  += c[mu]*cp[mu]*ss[mu]*ssp[mu]*A[mu][nu]*s2[nu];
      CxSPxSSPxAtxS2 += c[mu]*sp[mu]*ssp[mu]*A[mu][nu]*s2[nu];
      CxCPxAtxS2 += c[mu]*cp[mu]*A[mu][nu]*s2[nu];
      S2xSSxSSPxAtxS2  += s2[mu]*ss[mu]*ssp[mu]*A[mu][nu]*s2[nu];
      S2xSSP2xAtxS2  += s2[mu]*ssp2[mu]*A[mu][nu]*s2[nu];
      SP2xSSP2xAtxS2 += sp2[mu]*ssp2[mu]*A[mu][nu]*s2[nu];
      SSxSSPxAtxSSxSSP += ss[mu]*ssp[mu]*A[mu][nu]*ss[nu]*ssp[nu];
      SSxSSPxAtxSxSP += ss[mu]*ssp[mu]*A[mu][nu]*s[nu]*sp[nu];
      SSxSSPxAtxSS2  += ss[mu]*ssp[mu]*A[mu][nu]*ss2[nu];
      SSxSSPxAtxS2 += ss[mu]*ssp[mu]*A[mu][nu]*s2[nu];
      SxSPxAtxSxSP += s[mu]*sp[mu]*A[mu][nu]*s[nu]*sp[nu];
      SSP2xAtxS2 += ssp2[mu]*A[mu][nu]*s2[nu];
      SP2xAtxS2  += sp2[mu]*A[mu][nu]*s2[nu];

      PxCPxSxAtxCPxSxSSP += p[mu]*cp[mu]*s[mu]*A[mu][nu]*cp[nu]*s[nu]*ssp[nu];
      PxCPxSxAtxSxSP += p[mu]*cp[mu]*s[mu]*A[mu][nu]*s[nu]*sp[nu];
      PxCPxSPxAtxS2 += p[mu]*cp[mu]*sp[mu]*A[mu][nu]*s2[nu];
      PxCxCPxSSxAtxS2 += p[mu]*c[mu]*cp[mu]*ss[mu]*A[mu][nu]*s2[nu];
      PxCxSPxAtxS2 += p[mu]*c[mu]*sp[mu]*A[mu][nu]*s2[nu];
      PxSSPxAtxS2 += p[mu]*ssp[mu]*A[mu][nu]*s2[nu];
      PxSP2xSSPxAtxS2 += p[mu]*sp2[mu]*ssp[mu]*A[mu][nu]*s2[nu];
      PxS2xSSPxAtxS2 += p[mu]*s2[mu]*ssp[mu]*A[mu][nu]*s2[nu];
      PxS2xSSxAtxS2 += p[mu]*s2[mu]*ss[mu]*A[mu][nu]*s2[nu];
      PxSSxAtxCPxSxSS += p[mu]*ss[mu]*A[mu][nu]*cp[nu]*s[nu]*ss[nu];
      PxSSxAtxSSxSSP += p[mu]*ss[mu]*A[mu][nu]*ss[nu]*ssp[nu];
      PxSSxAtxSxSP += p[mu]*ss[mu]*A[mu][nu]*s[nu]*sp[nu];
      PxSSxAtxSS2 += p[mu]*ss[mu]*A[mu][nu]*ss2[nu];
      PxSSxAtxS2 += p[mu]*ss[mu]*A[mu][nu]*s2[nu];
    }
  }

  // asm("#CIAO!");

  Delta1 = S2;
  Delta1p = SQ2;
  Delta2 = SS2 + 4.0*S2*S2;
  Delta2p = SSP2 + 4.0*SQ2*SQ2;

  Delta12  = Delta1*Delta1;
  Delta13  = Delta1*Delta12;
  Delta14  = Delta12*Delta12;
  Delta22  = Delta2*Delta2;
  Delta1p2 = Delta1p*Delta1p;
  Delta2p2 = Delta2p*Delta2p;
  Delta12p2 = Delta1*Delta2p2;

  Den = (Delta12*Delta2p2);

  /* Vertex Integrals (Int(ainf)-Int(a0))*/
  auto c0 = Fac*((-CPxSxSSP*CPxSxSSP/(4.*Den) + SSP2/(4.*Delta12p2)) - (-((16.0*Delta12 - Delta22)*SS2*(-4.0*Delta1 + SS2))/(256.0*Delta14*Delta22)));
  Int[0] = c0;
  auto c1 = Fac*((CPxSxSSP*CPxSxSSP/(4.*Den))-(((-1.0/(16.*Delta12) + 1.0/(Delta22))*SS2*SS2)/(16.0*Delta12)));
  Int[1] = c1;
  auto c2 = Fac*((-CPxSxSSPxAtxCPxSxSSP/(4.*Den)+(CPxSPxSSPxAtxS2*Delta1p)/(Den)-(CPxSxSSPxAtxSxSP*Delta1p)/(Den)+(CPxSPxSSP*Delta1p)/(Delta12p2)+(Delta1p2*SP2)/(Delta12p2)+(Delta1p2*SP2xAtxS2)/(Den)-SP2xSSP2/(4.*Delta12p2)-SP2xSSP2xAtxS2/(4.*Den)+SSP2xAtxS2/(4.*Den)-(CPxSxSSP*Delta1p*SxSP)/(Den)-(Delta1p2*SxSP*SxSP)/(Den)-(Delta1p2*SxSPxAtxSxSP)/(Den))-(-(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.0*Delta12*Delta22)));
  Int[2] = c2;
  auto c3 = Fac*(((CPxSxSSP*Delta1p*SxSP)/(Den)+(Delta1p2*SxSP*SxSP)/(Den))-((2.0*S2*S2 + SS2)/(2.0*Delta22)));
  Int[3] = c3;
  auto c4 = Fac*(((CPxSPxSSPxAtxS2*Delta1p)/(4.*Den)-(CPxSxSSPxAtxSxSP*Delta1p)/(4.*Den)+(CPxSPxSSP*Delta1p)/(4.*Delta12p2)-(AtxS2*Delta1p2)/(4.*Den)-(3.0*Delta1p2)/(4.*Delta12p2)-(Delta1p2*S2xSP2)/(4.*Den)+(Delta1p2*SP2)/(4.*Delta12p2)+(Delta1p2*SP2xAtxS2)/(4.*Den)-(SP2*SSP2)/(16.*Delta12p2)-(SP2xAtxS2*SSP2)/(16.*Den)-(CPxSxSSP*Delta1p*SxSP)/(4.*Den)+(SSP2*SxSP*SxSP)/(16.*Den)+(SSP2*SxSPxAtxSxSP)/(16.*Den))-((-4.0*AtxS2 - 16.0*S2 + 4.0*S2*S2 + 4.0*S2xAtxS2 + SS2)/(16.0*Delta22)));
  Int[4] = c4;
  auto c5 = Fac*((-Delta1p2/(4.*Delta12p2) + (Delta1p2*S2xSP2)/(4.*Den) + (CPxSxSSP*Delta1p*SxSP)/(4.*Den) - (SSP2*SxSP*SxSP)/(16.*Den))-(0.0));
  Int[5] = c5;
  auto c6 = Fac*((CPxSxSSPxAtxCPxSxSSP/(4.*Den)+SP2xSSP2/(4.*Delta12p2)+SP2xSSP2xAtxS2/(4.*Den)+(AtxS2*SSP2)/(4.*Den)+(S2xSP2*SSP2)/(4.*Den)-(SP2*SSP2)/(4.*Delta12p2)-(SP2xAtxS2*SSP2)/(4.*Den)-SSP2xAtxS2/(4.*Den))-((4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 4.0*Delta1*SS2 - 4.0*Delta12*SS2 - 4.0*S2xAtxS2*SS2 - SS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.0*Delta12*Delta22)));
  Int[6] = c6;
  auto c7 = Fac*((-(S2xSP2*SSP2)/(4.*Den))-((SS2*(-4.0*Delta1 + SS2))/(16.0*Delta12*Delta22)));
  Int[7] = c7;
  auto c8 = Fac*(((CPxSxSSxAtxSSxSSP*Delta1p)/(16.*Den)-(CxCPxSSxSSPxAtxS2*Delta1p)/(8.*Den)-(CxCPxSSxSSP*Delta1p)/(8.*Delta12p2)+(CxCPxAtxS2*Delta1p*SSxSSP)/(8.*Den)+(CxCP*Delta1p*SSxSSP)/(8.*Delta12p2))-((4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 12.0*Delta1*SS2 - 4.0*Delta1*S2*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(32.0*Delta1*Delta22)));
  Int[8] = c8;
  auto c9 = Fac*(((S2xSSP2*SS2)/(64.*Delta12p2)+(S2xSSP2xAtxS2*SS2)/(64.*Den)-(SS2*SSP2)/(64.*Delta12p2)-(SS2*SSP2xAtxS2)/(64.*Den)-(S2xSSxSSP*SSxSSP)/(32.*Delta12p2)-(S2xSSxSSPxAtxS2*SSxSSP)/(32.*Den)+SSxSSP*SSxSSP/(64.*Delta2p2)-(AtxS2*SSxSSP*SSxSSP)/(64.*Den)-SSxSSP*SSxSSP/(32.*Delta12p2)+(S2xAtxS2*SSxSSP*SSxSSP)/(64.*Den)+(SSxSSP*SSxSSPxAtxS2)/(32.*Den)-(SSxSSP*SSxSSPxAtxSS2)/(128.*Den)+(SS2*SSxSSPxAtxSSxSSP)/(256.*Den))-(-(SS2*(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 12.0*Delta1*SS2 - 4.0*Delta12*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2))/(256.0*Delta12*Delta22)));
  Int[9] = c9;
  auto c10 = Fac*((-(S2xSSP2*SS2)/(16.*Delta12p2)-(S2xSSP2xAtxS2*SS2)/(16.*Den)-(S2xSS2*SSP2)/(16.*Delta12p2)-(S2xSS2xAtxS2*SSP2)/(16.*Den)+(SS2*SSP2)/(16.*Delta2p2)-(AtxS2*SS2*SSP2)/(16.*Den)-(SS2*SSP2)/(8.*Delta12p2)+(S2xAtxS2*SS2*SSP2)/(16.*Den)+(SS2xAtxS2*SSP2)/(16.*Den)-(SS2xAtxSS2*SSP2)/(64.*Den)+(SS2*SSP2xAtxS2)/(16.*Den)+(S2xSSxSSP*SSxSSP)/(8.*Delta12p2)+(S2xSSxSSPxAtxS2*SSxSSP)/(8.*Den)-SSxSSP*SSxSSP/(16.*Delta2p2)+(AtxS2*SSxSSP*SSxSSP)/(16.*Den)+SSxSSP*SSxSSP/(8.*Delta12p2)-(S2xAtxS2*SSxSSP*SSxSSP)/(16.*Den)-(SSxSSP*SSxSSPxAtxS2)/(8.*Den)+(SSxSSP*SSxSSPxAtxSS2)/(32.*Den)-(SS2*SSxSSPxAtxSSxSSP)/(64.*Den))-(0.0));
  Int[10] = c10;
  auto c11 = Fac*(((Delta1p2*S2xSS2)/(4.*Delta12p2)+(Delta1p2*S2xSS2xAtxS2)/(4.*Den)-(Delta1p2*SS2)/(4.*Delta2p2)+(AtxS2*Delta1p2*SS2)/(4.*Den)+(3.0*Delta1p2*SS2)/(4.*Delta12p2)-(Delta1p2*S2xAtxS2*SS2)/(4.*Den)-(Delta1p2*SS2xAtxS2)/(4.*Den)+(Delta1p2*SS2xAtxSS2)/(16.*Den))-((4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 - 4.0*S2*S2*SS2 - 4.0*S2xAtxS2*SS2 + 4.0*S2*(S2xSS2 + 3.0*SS2) - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.*Delta22)));
  Int[11] = c11;

  /* Self-Energy Integrals (IntS(ainf)-IntS(a0))*/
  auto d0 = Fac*((-(2.0*CPxSxSSP*PxCPxS + Delta1*PxSSP)/(4.0*Delta12*Delta2p*P2))-((Delta1*S2xSS2*((0.125*Delta1 - 0.03125*Delta2)*Delta2 - 0.125*Delta1*SS2) + SS2*(Delta1*Delta2*(-0.0625*Delta1 + 0.0625*Delta12 + 0.015625*Delta2 - 0.015625*Delta1*Delta2) + (0.0625*Delta12 + 0.125*Delta13 - 0.00390625*Delta22)*SS2))/(Delta14*Delta22)));
  IntS[0] = d0;
  auto d1 = Fac*(((2.0*CPxSxSSP*PxCPxS - Delta1*PxSSP)/(4.0*Delta12*Delta2p*P2))-((Delta1*S2xSS2*((-0.125*Delta1 + 0.03125*Delta2)*Delta2 + 0.125*Delta1*SS2) + SS2*(Delta12*(-0.0625*Delta1 + 0.015625*Delta2)*Delta2 + (-0.0625*Delta12 - 0.125*Delta13 + 0.00390625*Delta22)*SS2))/(Delta14*Delta22)));
  IntS[1] = d1;
  auto d2 = Fac*(((32.0*Delta1*Delta1p*PxCPxSP + 32.0*Delta1p*PxCPxSPxAtxS2 - 16.0*PxCPxSxAtxCPxSxSSP - 32.0*Delta1p*PxCPxSxAtxSxSP + AtxS2*Delta2p*PxSINP + 3.0*Delta1*Delta2p*PxSINP - 16.0*Delta1*PxSP2xSSP - 16.0*PxSP2xSSPxAtxS2 - 8.0*AtxS2*PxSSP + 16.0*PxSSPxAtxS2 - 8.0*PxSSP*S2xSP2 - 32.0*Delta1p*PxCPxS*SxSP + 8.0*PxSSP*SxSP*SxSP + 8.0*PxSSP*SxSPxAtxSxSP)/(32.*Delta12*Delta2p*P2))-((8.0*CCxS2xSS2xAtxS2 + 8.0*CCxS2xSS2*Delta1 - 8.0*Delta1*Delta2 + 48.0*Delta12*Delta2 + 3.0*Delta1*Delta22 - 12.0*Delta12*Delta2*S2 + 4.0*Delta2*S2xAtxS2 - 12.0*Delta1*Delta2*S2xAtxS2 - 8.0*Delta1*S2xSS2 + 24.0*Delta12*S2xSS2 - 2.0*Delta2*S2xSS2 + 8.0*S2xAtxS2*S2xSS2 + 16.0*S2xSS2xAtxS2 + 16.0*Delta1*S2xSS2xAtxS2 - 4.0*S2xSS2xAtxSS2 + 4.0*Delta1*SS2 + 4.0*Delta12*SS2 + 2.0*Delta2*SS2 - 7.0*Delta1*Delta2*SS2 - 8.0*Delta12*S2*SS2 - 4.0*S2xAtxS2*SS2 - 8.0*Delta1*S2xAtxS2*SS2 + 2.0*S2xSS2*SS2 - SS2*SS2 - 2.0*Delta1*SS2*SS2 + AtxS2*(Delta2*(-4.0 + 12.0*Delta1 + Delta2) - 8.0*S2xSS2 + (4.0 + 8.0*Delta1)*SS2)- 8.0*SS2xAtxS2 - 16.0*Delta1*SS2xAtxS2 - 2.0*Delta2*SS2xAtxS2 + 2.0*SS2xAtxSS2 + 4.0*Delta1*SS2xAtxSS2)/(32.*Delta12*Delta22)));
  IntS[2] = d2;
  auto d3 = Fac*(((Delta1*Delta2p*PxSINP + 8.0*PxSSP*S2xSP2 + 32.0*Delta1p*PxCPxS*SxSP - 8.0*PxSSP*SxSP*SxSP)/(32.*Delta12*Delta2p*P2))-((Delta1*Delta2*(8.0 + Delta2 - 4.0*Delta1*(3.0 + S2)) + S2xSS2*(8.0*Delta1 + 8.0*Delta12 + 2.0*Delta2 - 2.0*SS2) - (Delta1*(4.0 - 3.0*Delta2) + 2.0*Delta2 + 4.0*Delta12*(3.0 + 2.0*S2))*SS2 + (1.0 + 2.0*Delta1)*SS2*SS2)/(32.*Delta12*Delta22)));
  IntS[3] = d3;
  auto d4 = Fac*(((4.0*Delta1*Delta1p*PxCxCPxSS + 4.0*Delta1p*PxCxCPxSSxAtxS2 + 2.0*CxSPxSSPxAtxS2*PxSS + 2.0*CxSPxSSP*Delta1*PxSS - 4.0*CxCPxAtxS2*Delta1p*PxSS - 4.0*CxCP*Delta1*Delta1p*PxSS - 2.0*Delta1p*PxSSxAtxCPxSxSS - 2.0*Delta1*PxCxSP*SSxSSP - 2.0*PxCxSPxAtxS2*SSxSSP + PxSSxAtxSxSP*SSxSSP - PxSS*SSxSSPxAtxSxSP)/(8.*Delta12*Delta2p*P2))-((8.0*CCxS2xSS2xAtxS2*Delta1 + 8.0*CCxS2xSS2*Delta12-48.0*Delta12*S2xSS2-4.0*Delta1*Delta2*S2xSS2+32.0*Delta12*S2*S2xSS2+16.0*Delta1*S2xAtxS2*S2xSS2+16.0*Delta1*S2xSS2xAtxS2+16.0*Delta12*S2xSS2xAtxS2-4.0*Delta2*S2xSS2xAtxS2-4.0*Delta1*S2xSS2xAtxSS2+24.0*Delta12*SS2-12.0*Delta1*Delta2*SS2+4.0*Delta12*Delta2*SS2+40*Delta12*S2*SS2-16.0*Delta12*S2*S2*SS2-8.0*Delta1*S2xAtxS2*SS2-16.0*Delta12*S2xAtxS2*SS2+4.0*Delta2*S2xAtxS2*SS2-4.0*AtxS2*(4.0*Delta1*S2xSS2+(-2.0*Delta1-4.0*Delta12+Delta2)*SS2)-8.0*Delta1*SS2xAtxS2-16.0*Delta12*SS2xAtxS2+4.0*Delta2*SS2xAtxS2+2.0*Delta1*SS2xAtxSS2+4.0*Delta12*SS2xAtxSS2-Delta2*SS2xAtxSS2)/(32.*Delta12*Delta22)));
  IntS[4] = d4;
  auto d5 = Fac*(((-8.0*Delta1*PxS2xSSP*SS2 - 8.0*PxS2xSSPxAtxS2*SS2 + 8.0*PxSSPxAtxS2*SS2 - 2.0*PxSSxAtxSSxSSP*SS2 - PxSSP*(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 4.0*Delta1*SS2 - 4.0*Delta12*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2) + 8.0*Delta1*PxS2xSS*SSxSSP + 8.0*PxS2xSSxAtxS2*SSxSSP - 8.0*PxSSxAtxS2*SSxSSP + 2.0*PxSSxAtxSS2*SSxSSP + 2.0*PxSS*(4.0*Delta1*S2xSSxSSP + 4.0*S2xSSxSSPxAtxS2 + 4.0*AtxS2*SSxSSP + 8.0*Delta1*SSxSSP - 4.0*Delta12*SSxSSP - 4.0*S2xAtxS2*SSxSSP - 4.0*SSxSSPxAtxS2 + SSxSSPxAtxSS2))/(64.*Delta12*Delta2p*P2))-((8.0*CCxS2xSS2xAtxS2*Delta2 + 8.0*CCxS2xSS2*Delta1*Delta2+8.0*Delta1*S2xSS2*S2xSS2-4.0*Delta2*S2xSS2xAtxSS2+4.0*AtxS2*Delta1*Delta2*SS2+12.0*Delta12*Delta2*SS2-4.0*Delta12*Delta2*S2*SS2-4.0*Delta1*Delta2*S2xAtxS2*SS2-4.0*AtxS2*SS2*SS2-12.0*Delta1*SS2*SS2-8.0*AtxS2*Delta1*SS2*SS2-20.0*Delta12*SS2*SS2-2.0*Delta1*Delta2*SS2*SS2+8.0*Delta12*S2*SS2*SS2+4.0*S2xAtxS2*SS2*SS2+8.0*Delta1*S2xAtxS2*SS2*SS2-4.0*S2xSS2xAtxS2*(-((2.0+Delta1)*Delta2)+SS2+2.0*Delta1*SS2)-4.0*Delta1*Delta2*SS2xAtxS2+4.0*SS2*SS2xAtxS2+8.0*Delta1*SS2*SS2xAtxS2-2.0*Delta2*SS2*SS2xAtxS2+Delta1*Delta2*SS2xAtxSS2-SS2*SS2xAtxSS2-2.0*Delta1*SS2*SS2xAtxSS2+2.0*S2xSS2*(-4.0*AtxS2*Delta2-12.0*Delta1*Delta2+6.0*Delta12*Delta2+4.0*Delta2*S2xAtxS2+4.0*S2xSS2xAtxS2+4.0*AtxS2*SS2+10.0*Delta1*SS2-8.0*Delta12*SS2-4.0*S2xAtxS2*SS2-4.0*SS2xAtxS2+SS2xAtxSS2))/(128.*Delta12*Delta22)));
  IntS[5] = d5;
  // IntS[6]  = 0.0;
  // IntS[7]  = 0.0;
  // IntS[8]  = 0.0;
  // IntS[9]  = 0.0;

}


int main(int narg,char **arg)
{
  int nthreads;

  omp_set_max_active_levels(1);
  #pragma omp parallel
  #pragma omp master
  {
    system("clear");

    nthreads = omp_get_num_threads();
    cout<<"Using "<<nthreads<<" threads"<<endl<<endl;
  }

  char path_glb[128]="input.txt";

  if(narg>2){
    cerr<<"Number of arguments not valid."<<endl;
    exit(0);
  }
  if(narg==2){
    string path=arg[1];
    strcpy(path_glb,path.c_str());
  }

  /* Reading input file */
  cout<<"Reading input from \""<<path_glb<<"\"."<<endl;
  read_input_glb(path_glb);

  /* Loop momenta */
  int LT = L+T;
  j.resize(LT);
  kt2_v.resize(LT);
  k.resize(LT);
  kt.resize(LT);
  kt4_v.resize(LT);
  kt6_v.resize(LT);
  s_v.resize(LT);
  ss_v.resize(LT);
  s2_v.resize(LT);
  ss2_v.resize(LT);
  c_v.resize(LT);
  cc_v.resize(LT);

  cout<<"Vectorization with N = "<<N<<endl;
  cout<<endl;
  cout<<"Real type: "<<TO_STRING(REAL)<<endl;
  

  #pragma omp parallel for
  for(int l=0; l<LT; l++)
  {
    int p=(l>=L) /* 0 if l<L or 1 if l>=L */;
    int i=l-p*L;

    k[l] = 2.0*M_PI/dim[p]*(i+0.5) - al*sin(2.0*M_PI/dim[p]*(i+0.5));

    kt[l] = 2.0*sin(k[l]/2.0); /*ktilde*/
    kt2_v[l] = kt[l]*kt[l];
    kt4_v[l] = kt2_v[l]*kt2_v[l];
    kt6_v[l] = kt4_v[l]*kt2_v[l];

    j[l] = 1.0 - cos(2.0*M_PI/dim[p]*(i+0.5));

    s_v[l]    = sin(k[l]/2.0);
    ss_v[l]   = sin(k[l]);

    c_v[l]    = cos(k[l]/2.0);
    cc_v[l]   = cos(k[l]);

    s2_v[l]   = s_v[l]*s_v[l];
    ss2_v[l]  = ss_v[l]*ss_v[l];
  }

  find_eqmoms();

  // Compute RC corrections for independent external momenta
  vector<vector<Real>> DeltaZ_moms, DeltaG_moms;

  for(int imom=0; imom<eqmoms; imom++)
  {
    cout<<"\r"<<imom+1<<"/"<<eqmoms<<flush;

    array<Real,4> ap = ap_eq[imom];

    P2 = norm4(ap);
    PxSINP = P2;
    Np0 = 4.0 - count(ap.begin(), ap.end(), 0.0);
    
    Real eps=
#if SIMD_INST_SET == FLOAt128
      1.0e-28
#else
      1.0e-14
#endif
      ;
    for(int mu=0; mu<4; mu++)
    {
      if(fabs(sin(ap[mu]))>eps) p[mu] = P2/Np0/sin(ap[mu]);
      else p[mu] = 0.0;
    }

    int L3T = 3*L+T;
    kn.resize(L3T);
    kp.resize(L3T);
    sp_v.resize(L3T);
    cp_v.resize(L3T);
    sp2_v.resize(L3T);
    sq2_v.resize(L3T);
    ssp_v.resize(L3T);
    ssp2_v.resize(L3T);

    #pragma omp parallel for
    for(int j=0; j<L3T; j++)
    {
      int q = (j>=L)+(j>=2*L);
      int p = q+(j>=3*L);
      int l = j-q*L;

      kn[j] = k[l] + 2.0*ap[p];
      kp[j] = k[l] + ap[p];

      sp_v[j]   = sin(kn[j]/2.0);
      cp_v[j]   = cos(kn[j]/2.0);

      sp2_v[j]  = sp_v[j]*sp_v[j];
      sq2_v[j]  = pow(sin(kp[j]/2.0),2.0);

      ssp_v[j]  = sin(kp[j]);
      ssp2_v[j] = ssp_v[j]*ssp_v[j];
    }

    Real Int[12]={0.0};
    Real IntS[10]={0.0};

    if(L*L*L*T%N!=0)
    {
      cerr<<"Volume = "<<L*L*L*T<<" not divisible by SIMD size "<<N<<endl;
      exit(0);
    }

    #pragma omp parallel for reduction(+:Int,IntS)
    // #pragma omp parallel for
    for(int i=0;i<L*L*L*T/N;i++)
    {
      std::array<vtype,12> IntV;
      std::array<vtype,10> IntSV;

      std::array<int,N> u;
      for(int j=0;j<N;j++) u[j]=i*N+j;

      intComp(IntV,IntSV,u);

      for(int j=0;j<N;j++)
      {
        for(int o=0;o<12;o++) Int[o]+=IntV[o][j];
        for(int o=0;o<10;o++) IntS[o]+=IntSV[o][j];
      }
    }

    Int[0] -= (7.0 + M_PI*M_PI*(44.0*Z0 - 6.0))/6.0;
    Int[1] -= 4.0/3.0 - EulerGamma + F0 - log(P2) + M_PI*M_PI - 28.0*Z0*M_PI*M_PI/3.0;

    IntS[0] -= 2.0*(-1.0 + M_PI*M_PI*(-3.0+22.0*Z0))/3.0;
    IntS[1] -= -4.0/3.0 + EulerGamma - F0 + log(P2) + 2.0*M_PI*M_PI - 38.0*M_PI*M_PI*Z0/3.0;

    // printf("--- INTEGRALS -------------------- \n");
    // printf("i \t IntV \t \t IntV(a=0)\n");
    // printf("---------------------------------- \n");
    // for(int i=0; i<12; i++)
    // {
    //   printf("%d \t %lf",i,Int[i]);
    //   if(i<10) printf(" \t %lf\n",IntS[i]);
    //   else printf("\n");
    // }

    // string RCs[6] = {"q","S","P","V","A","T"};
    vector<Real> DeltaZ = compute_Z(Int,IntS);
    vector<Real> DeltaG = compute_Gamma(Int,IntS);


    // printf("--- CORRECTIONS O(g2ainf) -- \n");
    // printf("i \t DeltaZ\n");
    // printf("---------------------------- \n");
    // for(int i=0; i<6; i++)
    // {
    //   printf("%s \t %lf\n",RCs[i].c_str(),DeltaZ[i]);
    // }

    DeltaZ_moms.push_back(DeltaZ);
    DeltaG_moms.push_back(DeltaG);
  }

  vector<ofstream>  DeltaG_file(6), DeltaZ_file(6);
  string RCs[6] = {"q","S","P","V","A","T"};
  for(int iRC=0;iRC<6;iRC++)
  {
    DeltaG_file[iRC].open("DeltaG"+RCs[iRC]);
    DeltaZ_file[iRC].open("DeltaZ"+RCs[iRC]);

    int moms_max, ind_imom;
    if(eq) moms_max=eqmoms;  /* only eqmoms*/
    else   moms_max=moms;

    for(int imom=0; imom<moms_max; imom++)
    {
      if(eq) ind_imom = imom;
      else   ind_imom = tag_list[imom];
      
      DeltaG_file[iRC]<<(double)DeltaG_moms[ind_imom][iRC]<<endl;
      DeltaZ_file[iRC]<<(double)DeltaZ_moms[ind_imom][iRC]<<endl;
    }
  }

  exit(0);
}
