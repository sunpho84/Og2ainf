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
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "aliases.hpp"
#include "global.hpp"
#include "read_input.hpp"
#include "functions.hpp"

#define EulerGamma 0.5772156649015329
#define F0 4.36923
#define Z0 0.154933

#pragma omp declare reduction(vd_t_plus : vd_t : omp_out += omp_in) initializer(omp_priv=omp_orig)

Vector4d p;
double P2;
double PxSINP;
vector<double> j,k,kt,kt2_v,kt4_v,kt6_v,s_v,ss_v,s2_v,ss2_v,c_v,cc_v,
  kn,kp,sp_v,cp_v,sp2_v,sq2_v,ssp_v,ssp2_v;

void intComp(vd_t& Int, vd_t& Inta0,vd_t& IntS,vd_t& IntSa0,int i)
{
  int j=i;
  int i4=j%T;
  j/=T;
  int i3=j%L;
  j/=L;
  int i2=j%L;
  j/=L;
  int i1=j;

  Vector4d kt2, kt4, kt6;
  Vector4d s,ss,s2,ss2,c,cc;
  Vector4d sp,cp,sp2,sq2,ssp,ssp2;
  double S2, SP2, SQ2, SS2, SSP2;
  double Fac;
  double Delta1, Delta2, Delta1p, Delta2p;
  double Delta12, Delta13, Delta14, Delta22, Delta1p2, Delta2p2, Delta12p2, Den;
  double SxSP, CxCP, SSxSSP, S2xSP2, S2xSS2, S2xSSP2, SP2xSSP2, S2xSSxSSP, CPxSxSS, CPxSxSSP, CPxSPxSSP, CxCPxSSxSSP, CxSPxSSP, CCxS2xSS2;
  double PxSS, PxSSP, PxCPxS, PxS2xSSP, PxS2xSS, PxCxCPxSS, PxCPxSP, PxCxSP, PxSP2xSSP;
  double CCxS2xSS2xAtxS2, S2xSS2xAtxSS2, S2xSS2xAtxS2, SS2xAtxSS2, SS2xAtxS2, S2xAtxS2, AtxS2;
  double CPxSxSSPxAtxCPxSxSSP, CPxSxSSxAtxSSxSSP, CPxSxSSPxAtxSxSP, CPxSPxSSPxAtxS2, CxCPxSSxSSPxAtxS2, CxSPxSSPxAtxS2, CxCPxAtxS2, S2xSSxSSPxAtxS2, S2xSSP2xAtxS2, SP2xSSP2xAtxS2, SSxSSPxAtxSSxSSP, SSxSSPxAtxSxSP, SSxSSPxAtxSS2, SSxSSPxAtxS2, SxSPxAtxSxSP, SSP2xAtxS2, SP2xAtxS2;
  double PxCPxSxAtxCPxSxSSP, PxCPxSxAtxSxSP, PxCPxSPxAtxS2, PxCxCPxSSxAtxS2, PxCxSPxAtxS2, PxSSPxAtxS2, PxSP2xSSPxAtxS2, PxS2xSSPxAtxS2, PxS2xSSxAtxS2, PxSSxAtxCPxSxSS, PxSSxAtxSSxSSP, PxSSxAtxSxSP, PxSSxAtxSS2, PxSSxAtxS2;
  Matrix4d A;

  Fac = j[i1]*j[i2]*j[i3]*j[i4+L]*16.0*M_PI*M_PI/(L*L*L*T);

  kt2 << kt2_v[i1],kt2_v[i2],kt2_v[i3],kt2_v[i4+L];
  kt4 << kt4_v[i1],kt4_v[i2],kt4_v[i3],kt4_v[i4+L];
  kt6 << kt6_v[i1],kt6_v[i2],kt6_v[i3],kt6_v[i4+L];

  s   << s_v[i1],s_v[i2],s_v[i3],s_v[i4+L];
  ss  << ss_v[i1],ss_v[i2],ss_v[i3],ss_v[i4+L];
  s2  << s2_v[i1],s2_v[i2],s2_v[i3],s2_v[i4+L];
  ss2 << ss2_v[i1],ss2_v[i2],ss2_v[i3],ss2_v[i4+L];
  c   << c_v[i1],c_v[i2],c_v[i3],c_v[i4+L];
  cc  << cc_v[i1],cc_v[i2],cc_v[i3],cc_v[i4+L];

  sp   << sp_v[i1],sp_v[i2+L],sp_v[i3+2*L],sp_v[i4+3*L];
  cp   << cp_v[i1],cp_v[i2+L],cp_v[i3+2*L],cp_v[i4+3*L];
  sp2  << sp2_v[i1],sp2_v[i2+L],sp2_v[i3+2*L],sp2_v[i4+3*L];
  sq2  << sq2_v[i1],sq2_v[i2+L],sq2_v[i3+2*L],sq2_v[i4+3*L];
  ssp  << ssp_v[i1],ssp_v[i2+L],ssp_v[i3+2*L],ssp_v[i4+3*L];
  ssp2 << ssp2_v[i1],ssp2_v[i2+L],ssp2_v[i3+2*L],ssp2_v[i4+3*L];

  A = make_propagator(kt2,kt4,kt6);

  S2   = prod(s2);
  SP2  = prod(sp2);
  SQ2  = prod(sq2);
  SS2  = prod(ss2);
  SSP2 = prod(ssp2);

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

  SxSP = prod(s,sp);
  CxCP = prod(c,cp);
  SSxSSP = prod(ss,ssp);
  S2xSP2 = prod(s2,sp2);
  S2xSS2 = prod(s2,ss2);
  S2xSSP2 = prod(s2,ssp2);
  SP2xSSP2 = prod(sp2,ssp2);
  S2xSSxSSP = prod(s2,ss,ssp);
  CPxSxSS = prod(cp,s,ss);
  CPxSxSSP = prod(cp,s,ssp);
  CPxSPxSSP = prod(cp,sp,ssp);
  CxCPxSSxSSP = prod(c,cp,ss,ssp);
  CxSPxSSP = prod(c,sp,ssp);
  CCxS2xSS2 = prod(cc,s2,ss2);

  PxSS = prod(p,ss);
  PxSSP  = prod(p,ssp);
  PxCPxS = prod(p,cp,s);
  PxS2xSSP = prod(p,s2,ssp);
  PxS2xSS = prod(p,s2,ss);
  PxCxCPxSS = prod(p,c,cp,ss);
  PxCPxSP = prod(p,cp,sp);
  PxCxSP = prod(p,c,sp);
  PxSP2xSSP = prod(p,sp2,ssp);

  CCxS2xSS2xAtxS2 = prod(cc,s2,ss2,A*s2);
  S2xSS2xAtxSS2 = prod(s2,ss2,A*ss2);
  S2xSS2xAtxS2 = prod(s2,ss2,A*s2);
  SS2xAtxSS2 = prod(ss2,A*ss2);
  SS2xAtxS2 = prod(ss2,A*s2);
  S2xAtxS2 = prod(s2,A*s2);
  AtxS2 = prod(A*s2);

  CPxSxSSPxAtxCPxSxSSP = prod(cp,s,ssp,A*cwiseProduct(cp,s,ssp));
  CPxSxSSxAtxSSxSSP = prod(cp,s,ss,A*cwiseProduct(ss,ssp));
  CPxSxSSPxAtxSxSP = prod(cp,s,ssp,A*cwiseProduct(s,sp));
  CPxSPxSSPxAtxS2 = prod(cp,sp,ssp,A*s2);
  CxCPxSSxSSPxAtxS2 = prod(c,cp,ss,ssp,A*s2);
  CxSPxSSPxAtxS2 = prod(c,sp,ssp,A*s2);
  CxCPxAtxS2 = prod(c,cp,A*s2);
  S2xSSxSSPxAtxS2 = prod(s2,ss,ssp,A*s2);
  S2xSSP2xAtxS2 = prod(s2,ssp2,A*s2);
  SP2xSSP2xAtxS2 = prod(sp2,ssp2,A*s2);
  SSxSSPxAtxSSxSSP = prod(ss,ssp,A*cwiseProduct(ss,ssp));
  SSxSSPxAtxSxSP = prod(ss,ssp,A*cwiseProduct(s,sp));
  SSxSSPxAtxSS2 = prod(ss,ssp,A*ss2);
  SSxSSPxAtxS2 = prod(ss,ssp,A*s2);
  SxSPxAtxSxSP = prod(s,sp,A*cwiseProduct(s,sp));
  SSP2xAtxS2 = prod(ssp2,A*s2);
  SP2xAtxS2 = prod(sp2,A*s2);

  PxCPxSxAtxCPxSxSSP = prod(p,cp,s,A*cwiseProduct(cp,s,ssp));
  PxCPxSxAtxSxSP = prod(p,cp,s,A*cwiseProduct(s,sp));
  PxCPxSPxAtxS2 = prod(p,cp,sp,A*s2);
  PxCxCPxSSxAtxS2 = prod(p,c,cp,ss,A*s2);
  PxCxSPxAtxS2 = prod(p,c,sp,A*s2);
  PxSSPxAtxS2 = prod(p,ssp,A*s2);
  PxSP2xSSPxAtxS2 = prod(p,sp2,ssp,A*s2);
  PxS2xSSPxAtxS2 = prod(p,s2,ssp,A*s2);
  PxS2xSSxAtxS2 = prod(p,s2,ss,A*s2);
  PxSSxAtxCPxSxSS = prod(p,ss,A*cwiseProduct(cp,s,ss));
  PxSSxAtxSSxSSP = prod(p,ss,A*cwiseProduct(ss,ssp));
  PxSSxAtxSxSP = prod(p,ss,A*cwiseProduct(s,sp));
  PxSSxAtxSS2 = prod(p,ss,A*ss2);
  PxSSxAtxS2 = prod(p,ss,A*s2);

  /* Vertex Integrals */
  Int[0]  += Fac*(-CPxSxSSP*CPxSxSSP/(4.*Den) + SSP2/(4.*Delta12p2));
  Int[1]  += Fac*(CPxSxSSP*CPxSxSSP/(4.*Den));
  Int[2]  += Fac*(-CPxSxSSPxAtxCPxSxSSP/(4.*Den)+(CPxSPxSSPxAtxS2*Delta1p)/(Den)-(CPxSxSSPxAtxSxSP*Delta1p)/(Den)+(CPxSPxSSP*Delta1p)/(Delta12p2)+(Delta1p2*SP2)/(Delta12p2)+(Delta1p2*SP2xAtxS2)/(Den)-SP2xSSP2/(4.*Delta12p2)-SP2xSSP2xAtxS2/(4.*Den)+SSP2xAtxS2/(4.*Den)-(CPxSxSSP*Delta1p*SxSP)/(Den)-(Delta1p2*SxSP*SxSP)/(Den)-(Delta1p2*SxSPxAtxSxSP)/(Den));
  Int[3]  += Fac*((CPxSxSSP*Delta1p*SxSP)/(Den)+(Delta1p2*SxSP*SxSP)/(Den));
  Int[4]  += Fac*((CPxSPxSSPxAtxS2*Delta1p)/(4.*Den)-(CPxSxSSPxAtxSxSP*Delta1p)/(4.*Den)+(CPxSPxSSP*Delta1p)/(4.*Delta12p2)-(AtxS2*Delta1p2)/(4.*Den)-(3*Delta1p2)/(4.*Delta12p2)-(Delta1p2*S2xSP2)/(4.*Den)+(Delta1p2*SP2)/(4.*Delta12p2)+(Delta1p2*SP2xAtxS2)/(4.*Den)-(SP2*SSP2)/(16.*Delta12p2)-(SP2xAtxS2*SSP2)/(16.*Den)-(CPxSxSSP*Delta1p*SxSP)/(4.*Den)+(SSP2*SxSP*SxSP)/(16.*Den)+(SSP2*SxSPxAtxSxSP)/(16.*Den));
  Int[5]  += Fac*(-Delta1p2/(4.*Delta12p2) + (Delta1p2*S2xSP2)/(4.*Den) + (CPxSxSSP*Delta1p*SxSP)/(4.*Den) - (SSP2*SxSP*SxSP)/(16.*Den));
  Int[6]  += Fac*(CPxSxSSPxAtxCPxSxSSP/(4.*Den)+SP2xSSP2/(4.*Delta12p2)+SP2xSSP2xAtxS2/(4.*Den)+(AtxS2*SSP2)/(4.*Den)+(S2xSP2*SSP2)/(4.*Den)-(SP2*SSP2)/(4.*Delta12p2)-(SP2xAtxS2*SSP2)/(4.*Den)-SSP2xAtxS2/(4.*Den));
  Int[7]  += Fac*(-(S2xSP2*SSP2)/(4.*Den));
  Int[8]  += Fac*((CPxSxSSxAtxSSxSSP*Delta1p)/(16.*Den)-(CxCPxSSxSSPxAtxS2*Delta1p)/(8.*Den)-(CxCPxSSxSSP*Delta1p)/(8.*Delta12p2)+(CxCPxAtxS2*Delta1p*SSxSSP)/(8.*Den)+(CxCP*Delta1p*SSxSSP)/(8.*Delta12p2));
  Int[9]  += Fac*((S2xSSP2*SS2)/(64.*Delta12p2)+(S2xSSP2xAtxS2*SS2)/(64.*Den)-(SS2*SSP2)/(64.*Delta12p2)-(SS2*SSP2xAtxS2)/(64.*Den)-(S2xSSxSSP*SSxSSP)/(32.*Delta12p2)-(S2xSSxSSPxAtxS2*SSxSSP)/(32.*Den)+SSxSSP*SSxSSP/(64.*Delta2p2)-(AtxS2*SSxSSP*SSxSSP)/(64.*Den)-SSxSSP*SSxSSP/(32.*Delta12p2)+(S2xAtxS2*SSxSSP*SSxSSP)/(64.*Den)+(SSxSSP*SSxSSPxAtxS2)/(32.*Den)-(SSxSSP*SSxSSPxAtxSS2)/(128.*Den)+(SS2*SSxSSPxAtxSSxSSP)/(256.*Den));
  Int[10] += Fac*(-(S2xSSP2*SS2)/(16.*Delta12p2)-(S2xSSP2xAtxS2*SS2)/(16.*Den)-(S2xSS2*SSP2)/(16.*Delta12p2)-(S2xSS2xAtxS2*SSP2)/(16.*Den)+(SS2*SSP2)/(16.*Delta2p2)-(AtxS2*SS2*SSP2)/(16.*Den)-(SS2*SSP2)/(8.*Delta12p2)+(S2xAtxS2*SS2*SSP2)/(16.*Den)+(SS2xAtxS2*SSP2)/(16.*Den)-(SS2xAtxSS2*SSP2)/(64.*Den)+(SS2*SSP2xAtxS2)/(16.*Den)+(S2xSSxSSP*SSxSSP)/(8.*Delta12p2)+(S2xSSxSSPxAtxS2*SSxSSP)/(8.*Den)-SSxSSP*SSxSSP/(16.*Delta2p2)+(AtxS2*SSxSSP*SSxSSP)/(16.*Den)+SSxSSP*SSxSSP/(8.*Delta12p2)-(S2xAtxS2*SSxSSP*SSxSSP)/(16.*Den)-(SSxSSP*SSxSSPxAtxS2)/(8.*Den)+(SSxSSP*SSxSSPxAtxSS2)/(32.*Den)-(SS2*SSxSSPxAtxSSxSSP)/(64.*Den));
  Int[11] += Fac*((Delta1p2*S2xSS2)/(4.*Delta12p2)+(Delta1p2*S2xSS2xAtxS2)/(4.*Den)-(Delta1p2*SS2)/(4.*Delta2p2)+(AtxS2*Delta1p2*SS2)/(4.*Den)+(3*Delta1p2*SS2)/(4.*Delta12p2)-(Delta1p2*S2xAtxS2*SS2)/(4.*Den)-(Delta1p2*SS2xAtxS2)/(4.*Den)+(Delta1p2*SS2xAtxSS2)/(16.*Den));

  /* Vertex Integrals a->0 */
  Inta0[0]  += Fac*(-((16.0*Delta12 - Delta22)*SS2*(-4.0*Delta1 + SS2))/(256.0*Delta14*Delta22));
  Inta0[1]  += Fac*(((-1.0/(16.*Delta12) + 1.0/(Delta22))*SS2*SS2)/(16.0*Delta12));
  Inta0[2]  += Fac*(-(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.0*Delta12*Delta22));
  Inta0[3]  += Fac*((2.0*S2*S2 + SS2)/(2.0*Delta22));
  Inta0[4]  += Fac*((-4.0*AtxS2 - 16.0*S2 + 4.0*S2*S2 + 4.0*S2xAtxS2 + SS2)/(16.0*Delta22));
  // Inta0[5]  += 0.0;
  Inta0[6]  += Fac*((4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 4.0*Delta1*SS2 - 4.0*Delta12*SS2 - 4.0*S2xAtxS2*SS2 - SS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.0*Delta12*Delta22));
  Inta0[7]  += Fac*((SS2*(-4.0*Delta1 + SS2))/(16.0*Delta12*Delta22));
  Inta0[8]  += Fac*((4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 12.0*Delta1*SS2 - 4.0*Delta1*S2*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(32.0*Delta1*Delta22));
  Inta0[9]  += Fac*(-(SS2*(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 12*Delta1*SS2 - 4.0*Delta12*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2))/(256.0*Delta12*Delta22));
  // Inta0[10] += 0.0;
  Inta0[11] += Fac*((4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 - 4.0*S2*S2*SS2 - 4.0*S2xAtxS2*SS2 + 4.0*S2*(S2xSS2 + 3*SS2) - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.*Delta22));

  IntS[0]  += Fac*(-(2.0*CPxSxSSP*PxCPxS + Delta1*PxSSP)/(4.0*Delta12*Delta2p*P2));
  IntS[1]  += Fac*((2.0*CPxSxSSP*PxCPxS - Delta1*PxSSP)/(4.0*Delta12*Delta2p*P2));
  IntS[2]  += Fac*((32.0*Delta1*Delta1p*PxCPxSP + 32*Delta1p*PxCPxSPxAtxS2 - 16.0*PxCPxSxAtxCPxSxSSP - 32*Delta1p*PxCPxSxAtxSxSP + AtxS2*Delta2p*PxSINP + 3*Delta1*Delta2p*PxSINP - 16.0*Delta1*PxSP2xSSP - 16.0*PxSP2xSSPxAtxS2 - 8.0*AtxS2*PxSSP + 16.0*PxSSPxAtxS2 - 8.0*PxSSP*S2xSP2 - 32*Delta1p*PxCPxS*SxSP + 8.0*PxSSP*SxSP*SxSP + 8.0*PxSSP*SxSPxAtxSxSP)/(32.*Delta12*Delta2p*P2));
  IntS[3]  += Fac*((Delta1*Delta2p*PxSINP + 8.0*PxSSP*S2xSP2 + 32*Delta1p*PxCPxS*SxSP - 8.0*PxSSP*SxSP*SxSP)/(32.*Delta12*Delta2p*P2));
  IntS[4]  += Fac*((4.0*Delta1*Delta1p*PxCxCPxSS + 4.0*Delta1p*PxCxCPxSSxAtxS2 + 2*CxSPxSSPxAtxS2*PxSS + 2*CxSPxSSP*Delta1*PxSS - 4.0*CxCPxAtxS2*Delta1p*PxSS - 4.0*CxCP*Delta1*Delta1p*PxSS - 2*Delta1p*PxSSxAtxCPxSxSS - 2*Delta1*PxCxSP*SSxSSP - 2*PxCxSPxAtxS2*SSxSSP + PxSSxAtxSxSP*SSxSSP - PxSS*SSxSSPxAtxSxSP)/(8.*Delta12*Delta2p*P2));
  IntS[5]  += Fac*((-8.0*Delta1*PxS2xSSP*SS2 - 8.0*PxS2xSSPxAtxS2*SS2 + 8.0*PxSSPxAtxS2*SS2 - 2.0*PxSSxAtxSSxSSP*SS2 - PxSSP*(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 4.0*Delta1*SS2 - 4.0*Delta12*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2) + 8.0*Delta1*PxS2xSS*SSxSSP + 8.0*PxS2xSSxAtxS2*SSxSSP - 8.0*PxSSxAtxS2*SSxSSP + 2*PxSSxAtxSS2*SSxSSP + 2*PxSS*(4.0*Delta1*S2xSSxSSP + 4.0*S2xSSxSSPxAtxS2 + 4.0*AtxS2*SSxSSP + 8.0*Delta1*SSxSSP - 4.0*Delta12*SSxSSP - 4.0*S2xAtxS2*SSxSSP - 4.0*SSxSSPxAtxS2 + SSxSSPxAtxSS2))/(64.*Delta12*Delta2p*P2));
  // IntS[6]  += 0.0;
  // IntS[7]  += 0.0;
  // IntS[8]  += 0.0;
  // IntS[9]  += 0.0;

  IntSa0[0]  += Fac*((Delta1*S2xSS2*((0.125*Delta1 - 0.03125*Delta2)*Delta2 - 0.125*Delta1*SS2) + SS2*(Delta1*Delta2*(-0.0625*Delta1 + 0.0625*Delta12 + 0.015625*Delta2 - 0.015625*Delta1*Delta2) + (0.0625*Delta12 + 0.125*Delta13 - 0.00390625*Delta22)*SS2))/(Delta14*Delta22));
  IntSa0[1]  += Fac*((Delta1*S2xSS2*((-0.125*Delta1 + 0.03125*Delta2)*Delta2 + 0.125*Delta1*SS2) + SS2*(Delta12*(-0.0625*Delta1 + 0.015625*Delta2)*Delta2 + (-0.0625*Delta12 - 0.125*Delta13 + 0.00390625*Delta22)*SS2))/(Delta14*Delta22));
  IntSa0[2]  += Fac*((8.0*CCxS2xSS2xAtxS2 + 8.0*CCxS2xSS2*Delta1 - 8.0*Delta1*Delta2 + 48.0*Delta12*Delta2 + 3*Delta1*Delta22 - 12*Delta12*Delta2*S2 + 4.0*Delta2*S2xAtxS2 - 12*Delta1*Delta2*S2xAtxS2 - 8.0*Delta1*S2xSS2 + 24.0*Delta12*S2xSS2 - 2*Delta2*S2xSS2 + 8.0*S2xAtxS2*S2xSS2 + 16.0*S2xSS2xAtxS2 + 16.0*Delta1*S2xSS2xAtxS2 - 4.0*S2xSS2xAtxSS2 + 4.0*Delta1*SS2 + 4.0*Delta12*SS2 + 2*Delta2*SS2 - 7*Delta1*Delta2*SS2 - 8.0*Delta12*S2*SS2 - 4.0*S2xAtxS2*SS2 - 8.0*Delta1*S2xAtxS2*SS2 + 2*S2xSS2*SS2 - SS2*SS2 - 2*Delta1*SS2*SS2 + AtxS2*(Delta2*(-4.0 + 12*Delta1 + Delta2) - 8.0*S2xSS2 + (4.0 + 8.0*Delta1)*SS2)- 8.0*SS2xAtxS2 - 16.0*Delta1*SS2xAtxS2 - 2*Delta2*SS2xAtxS2 + 2*SS2xAtxSS2 + 4.0*Delta1*SS2xAtxSS2)/(32.*Delta12*Delta22));
  IntSa0[3]  += Fac*((Delta1*Delta2*(8 + Delta2 - 4.0*Delta1*(3 + S2)) + S2xSS2*(8.0*Delta1 + 8.0*Delta12 + 2.0*Delta2 - 2.0*SS2) - (Delta1*(4.0 - 3.0*Delta2) + 2.0*Delta2 + 4.0*Delta12*(3 + 2*S2))*SS2 + (1 + 2*Delta1)*SS2*SS2)/(32.*Delta12*Delta22));
  IntSa0[4]  += Fac*((8.0*CCxS2xSS2xAtxS2*Delta1 + 8.0*CCxS2xSS2*Delta12-48.0*Delta12*S2xSS2-4.0*Delta1*Delta2*S2xSS2+32*Delta12*S2*S2xSS2+16.0*Delta1*S2xAtxS2*S2xSS2+16.0*Delta1*S2xSS2xAtxS2+16.0*Delta12*S2xSS2xAtxS2-4.0*Delta2*S2xSS2xAtxS2-4.0*Delta1*S2xSS2xAtxSS2+24.0*Delta12*SS2-12*Delta1*Delta2*SS2+4.0*Delta12*Delta2*SS2+40*Delta12*S2*SS2-16.0*Delta12*S2*S2*SS2-8.0*Delta1*S2xAtxS2*SS2-16.0*Delta12*S2xAtxS2*SS2+4.0*Delta2*S2xAtxS2*SS2-4.0*AtxS2*(4.0*Delta1*S2xSS2+(-2*Delta1-4.0*Delta12+Delta2)*SS2)-8.0*Delta1*SS2xAtxS2-16.0*Delta12*SS2xAtxS2+4.0*Delta2*SS2xAtxS2+2*Delta1*SS2xAtxSS2+4.0*Delta12*SS2xAtxSS2-Delta2*SS2xAtxSS2)/(32.*Delta12*Delta22));
  IntSa0[5]  += Fac*((8.0*CCxS2xSS2xAtxS2*Delta2 + 8.0*CCxS2xSS2*Delta1*Delta2+8.0*Delta1*S2xSS2*S2xSS2-4.0*Delta2*S2xSS2xAtxSS2+4.0*AtxS2*Delta1*Delta2*SS2+12*Delta12*Delta2*SS2-4.0*Delta12*Delta2*S2*SS2-4.0*Delta1*Delta2*S2xAtxS2*SS2-4.0*AtxS2*SS2*SS2-12*Delta1*SS2*SS2-8.0*AtxS2*Delta1*SS2*SS2-20*Delta12*SS2*SS2-2*Delta1*Delta2*SS2*SS2+8.0*Delta12*S2*SS2*SS2+4.0*S2xAtxS2*SS2*SS2+8.0*Delta1*S2xAtxS2*SS2*SS2-4.0*S2xSS2xAtxS2*(-((2.0+Delta1)*Delta2)+SS2+2.0*Delta1*SS2)-4.0*Delta1*Delta2*SS2xAtxS2+4.0*SS2*SS2xAtxS2+8.0*Delta1*SS2*SS2xAtxS2-2*Delta2*SS2*SS2xAtxS2+Delta1*Delta2*SS2xAtxSS2-SS2*SS2xAtxSS2-2*Delta1*SS2*SS2xAtxSS2+2*S2xSS2*(-4.0*AtxS2*Delta2-12*Delta1*Delta2+6.0*Delta12*Delta2+4.0*Delta2*S2xAtxS2+4.0*S2xSS2xAtxS2+4.0*AtxS2*SS2+10*Delta1*SS2-8.0*Delta12*SS2-4.0*S2xAtxS2*SS2-4.0*SS2xAtxS2+SS2xAtxSS2))/(128.*Delta12*Delta22));
  // IntSa0[6]  += 0.0;
  // IntSa0[7]  += 0.0;
  // IntSa0[8]  += 0.0;
  // IntSa0[9]  += 0.0;
}

int main(int narg,char **arg)
{
  omp_set_nested(1);
  #pragma omp parallel
  #pragma omp master
  {
    system("clear");

    cout<<"Using "<<omp_get_num_threads()<<" threads"<<endl<<endl;
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

  #pragma omp parallel for
  for(int p=0; p<dim.size(); p++){
    for(int i=0; i<dim(p); i++){

      int l = i + dim(0)*p;

      k[l] = 2.0*M_PI/dim(p)*(i+0.5) - al*sin(2.0*M_PI/dim(p)*(i+0.5));

      kt[l] = 2.0*sin(k[l]/2.0); /*ktilde*/
      kt2_v[l] = kt[l]*kt[l];
      kt4_v[l] = kt2_v[l]*kt2_v[l];
      kt6_v[l] = kt4_v[l]*kt2_v[l];

      j[l] = 1.0 - cos(2.0*M_PI/dim(p)*(i+0.5));

      s_v[l]    = sin(k[l]/2.0);
      ss_v[l]   = sin(k[l]);

      c_v[l]    = cos(k[l]/2.0);
      cc_v[l]   = cos(k[l]);

      s2_v[l]   = s_v[l]*s_v[l];
      ss2_v[l]  = ss_v[l]*ss_v[l];
    }
  }

  /* External momentum   --->  va in un loop */
  Vector4d np;
  np(0) = 0.0;
  np(1) = 4.0;
  np(2) = 5.0;
  np(3) = 6.0;
  if(APBC) np(3) += 0.5;

  Np0 = 4.0 - (np.array() == 0).count();

  Vector4d ap = 2.0*M_PI*np.cwiseProduct(V.cwiseInverse());
  double P2 = prod(ap,ap);

  Vector4d p = Vector4d::Zero();
  for(int i=0; i<4; i++)
  {
    if(sin(ap(i))) p(i) = P2/Np0/sin(ap(i));
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
  for(int p=0; p<V.size(); p++){
    for(int i=0; i<V(p); i++){
      int j = i + L*p;

      kn[j] = k[i + L*(int)(p/3)] + 2.0*ap(p);
      kp[j] = k[i + L*(int)(p/3)] + ap(p);

      sp_v[j]   = sin(kn[j]/2.0);
      cp_v[j]   = cos(kn[j]/2.0);

      sp2_v[j]  = sp_v[j]*sp_v[j];
      sq2_v[j]  = pow(sin(kp[j]/2.0),2.0);

      ssp_v[j]  = sin(kp[j]);
      ssp2_v[j] = ssp_v[j]*ssp_v[j];
    }
  }

 PxSINP = P2;

  vd_t Int(0.0,12), Inta0(0.0,12), IntS(0.0,10), IntSa0(0.0,10);

  #pragma omp parallel  for reduction(vd_t_plus:Int,Inta0,IntSa0,IntS)
  for(int i=0;i<L*L*L*T;i++)
    {
      intComp(Int,Inta0,IntS,IntSa0,i);
    }

  Inta0[0] += (7.0 + M_PI*M_PI*(44.0*Z0 - 6.0))/6.0;
  Inta0[1] += 4.0/3.0 - EulerGamma + F0 - log(P2) + M_PI*M_PI - 28.0*Z0*M_PI*M_PI/3.0;

  IntSa0[0] += 2.0*(-1.0 + M_PI*M_PI*(-3.0+22.0*Z0))/3.0;
  IntSa0[1] += -4.0/3.0 + EulerGamma - F0 + log(P2) + 2.0*M_PI*M_PI - 38.0*M_PI*M_PI*Z0/3.0;

  printf("--- INTEGRALS ---------------------------------------------------- \n");
  printf("i \t IntV \t \t IntV(a=0) \t IntS \t \t IntS(a=0)\n");
  printf("------------------------------------------------------------------ \n");
  for(int i=0; i<12; i++)
    {
      printf("%d \t %lf \t %lf",i,Int[i],Inta0[i]);
      if(i<10) printf(" \t %lf \t %lf\n",IntS[i],IntSa0[i]);
      else printf("\n");
    }

  // vector<double> Vert   = compute_V(Int);
  // vector<double> Verta0 = compute_V(Inta0);
  // double Sigma1   = compute_S(IntS);
  // double Sigma1a0 = compute_S(IntSa0);




  exit(0);
}
