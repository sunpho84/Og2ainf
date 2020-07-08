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
  vector<double> k(LT), kt(LT), j(LT);
  vector<double> kt2_v(LT), kt4_v(LT), kt6_v(LT);
  vector<double> s_v(LT), ss_v(LT), s2_v(LT), ss2_v(LT), c_v(LT), cc_v(LT);

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

      s2_v[l]   = pow(sin(k[l]/2.0),2.0);
      ss2_v[l]  = pow(sin(k[l]),2.0);
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
  Vector4d p = Vector4d::Zero();
  for(int i=0; i<4; i++)
  {
    if(sin(ap(i))) p(i) = 1.0/sin(ap(i));
  }

  int L3T = 3*L+T;
  vector<double> kn(L3T), kp(L3T);
  vector<double> sp_v(L3T), cp_v(L3T), sp2_v(L3T), sq2_v(L3T), ssp_v(L3T), ssp2_v(L3T);

  #pragma omp parallel for
  for(int p=0; p<V.size(); p++){
    for(int i=0; i<V(p); i++){
      int j = i + L*p;

      kn[j] = k[i + L*(int)(p/3)] + 2.0*ap(p);
      kp[j] = k[i + L*(int)(p/3)] + ap(p);

      sp_v[j]   = sin(kn[j]/2.0);
      cp_v[j]   = cos(kn[j]/2.0);

      sp2_v[j]  = pow(sin(kn[j]/2.0),2.0);
      sq2_v[j]  = pow(sin(kp[j]/2.0),2.0);

      ssp_v[j]  = sin(kp[j]);
      ssp2_v[j] = pow(sin(kp[j]),2.0);
    }
  }

  Vector4d kt2, kt4, kt6;
  Vector4d s,ss,s2,ss2,c,cc;
  Vector4d sp,cp,sp2,sq2,ssp,ssp2;
  double P2, S2, SP2, SQ2, SS2, SSP2;
  double Jac, Fac;
  double Delta1, Delta2, Delta1p, Delta2p, Den;
  double SxSP, CxCP, SSxSSP, S2xSP2, S2xSS2, S2xSSP2, SP2xSSP2, S2xSSxSSP, CPxSxSS, CPxSxSSP, CPxSPxSSP, CxCPxSSxSSP, CxSPxSSP, CCxS2xSS2;
  double PxSINP, PxSS, PxSSP, PxCPxS, PxS2xSSP, PxS2xSS, PxCxCPxSS, PxCPxSP, PxCxSP, PxSP2xSSP;
  double CCxS2xSS2xAtxS2, S2xSS2xAtxSS2, S2xSS2xAtxS2, SS2xAtxSS2, SS2xAtxS2, S2xAtxS2, AtxS2;
  double CPxSxSSPxAtxCPxSxSSP, CPxSxSSxAtxSSxSSP, CPxSxSSPxAtxSxSP, CPxSPxSSPxAtxS2, CxCPxSSxSSPxAtxS2, CxSPxSSPxAtxS2, CxCPxAtxS2, S2xSSxSSPxAtxS2, S2xSSP2xAtxS2, SP2xSSP2xAtxS2, SSxSSPxAtxSSxSSP, SSxSSPxAtxSxSP, SSxSSPxAtxSS2, SSxSSPxAtxS2, SxSPxAtxSxSP, SSP2xAtxS2, SP2xAtxS2;
  double PxCPxSxAtxCPxSxSSP, PxCPxSxAtxSxSP, PxCPxSPxAtxS2, PxCxCPxSSxAtxS2, PxCxSPxAtxS2, PxSSPxAtxS2, PxSP2xSSPxAtxS2, PxS2xSSPxAtxS2, PxS2xSSxAtxS2, PxSSxAtxCPxSxSS, PxSSxAtxSSxSSP, PxSSxAtxSxSP, PxSSxAtxSS2, PxSSxAtxS2;
  Matrix4d A;

  P2 = PxSINP = prod(ap,ap);

  vd_t Int(0.0,12), Inta0(0.0,12), IntS(0.0,10), IntSa0(0.0,10);

  Inta0[0] = (7.0 + M_PI*M_PI*(44.0*Z0 - 6.0))/6.0;
  Inta0[1] = 4.0/3.0 - EulerGamma + F0 - log(P2) + M_PI*M_PI - 28.0*Z0*M_PI*M_PI/3.0;

  IntSa0[0] = 2.0*(-1.0 + M_PI*M_PI*(-3.0+22.0*Z0))/3.0;
  IntSa0[1] = -4.0/3.0 + EulerGamma - F0 + log(P2) + 2.0*M_PI*M_PI - 38.0*M_PI*M_PI*Z0/3.0;

#pragma omp parallel for collapse(4)
  for(int i1=0; i1<L; i1++)
  for(int i2=0; i2<L; i2++)
  for(int i3=0; i3<L; i3++)
  for(int i4=0; i4<T; i4++)
  {
    Jac = j[i1]*j[i2]*j[i3]*j[14+L];

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
    Den = (Delta1*Delta1*Delta2p*Delta2p);

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

    Fac = Jac*16.0*M_PI*M_PI/(L*L*L*T);

    /* Vertex Integrals */
    Int[0]  += Fac*(-CPxSxSSP*CPxSxSSP/(4.*Den) + SSP2/(4.*Delta1*Delta2p*Delta2p));
    Int[1]  += Fac*(CPxSxSSP*CPxSxSSP/(4.*Den));
    Int[2]  += Fac*(-CPxSxSSPxAtxCPxSxSSP/(4.*Den)+(CPxSPxSSPxAtxS2*Delta1p)/(Den)-(CPxSxSSPxAtxSxSP*Delta1p)/(Den)+(CPxSPxSSP*Delta1p)/(Delta1*Delta2p*Delta2p)+(Delta1p*Delta1p*SP2)/(Delta1*Delta2p*Delta2p)+(Delta1p*Delta1p*SP2xAtxS2)/(Den)-SP2xSSP2/(4.*Delta1*Delta2p*Delta2p)-SP2xSSP2xAtxS2/(4.*Den)+SSP2xAtxS2/(4.*Den)-(CPxSxSSP*Delta1p*SxSP)/(Den)-(Delta1p*Delta1p*SxSP*SxSP)/(Den)-(Delta1p*Delta1p*SxSPxAtxSxSP)/(Den));
    Int[3]  += Fac*((CPxSxSSP*Delta1p*SxSP)/(Den)+(Delta1p*Delta1p*SxSP*SxSP)/(Den));
    Int[4]  += Fac*((CPxSPxSSPxAtxS2*Delta1p)/(4.*Den)-(CPxSxSSPxAtxSxSP*Delta1p)/(4.*Den)+(CPxSPxSSP*Delta1p)/(4.*Delta1*Delta2p*Delta2p)-(AtxS2*Delta1p*Delta1p)/(4.*Den)-(3*Delta1p*Delta1p)/(4.*Delta1*Delta2p*Delta2p)-(Delta1p*Delta1p*S2xSP2)/(4.*Den)+(Delta1p*Delta1p*SP2)/(4.*Delta1*Delta2p*Delta2p)+(Delta1p*Delta1p*SP2xAtxS2)/(4.*Den)-(SP2*SSP2)/(16.*Delta1*Delta2p*Delta2p)-(SP2xAtxS2*SSP2)/(16.*Den)-(CPxSxSSP*Delta1p*SxSP)/(4.*Den)+(SSP2*SxSP*SxSP)/(16.*Den)+(SSP2*SxSPxAtxSxSP)/(16.*Den));
    Int[5]  += Fac*(-Delta1p*Delta1p/(4.*Delta1*Delta2p*Delta2p) + (Delta1p*Delta1p*S2xSP2)/(4.*Den) + (CPxSxSSP*Delta1p*SxSP)/(4.*Den) - (SSP2*SxSP*SxSP)/(16.*Den));
    Int[6]  += Fac*(CPxSxSSPxAtxCPxSxSSP/(4.*Den)+SP2xSSP2/(4.*Delta1*Delta2p*Delta2p)+SP2xSSP2xAtxS2/(4.*Den)+(AtxS2*SSP2)/(4.*Den)+(S2xSP2*SSP2)/(4.*Den)-(SP2*SSP2)/(4.*Delta1*Delta2p*Delta2p)-(SP2xAtxS2*SSP2)/(4.*Den)-SSP2xAtxS2/(4.*Den));
    Int[7]  += Fac*(-(S2xSP2*SSP2)/(4.*Den));
    Int[8]  += Fac*((CPxSxSSxAtxSSxSSP*Delta1p)/(16.*Den)-(CxCPxSSxSSPxAtxS2*Delta1p)/(8.*Den)-(CxCPxSSxSSP*Delta1p)/(8.*Delta1*Delta2p*Delta2p)+(CxCPxAtxS2*Delta1p*SSxSSP)/(8.*Den)+(CxCP*Delta1p*SSxSSP)/(8.*Delta1*Delta2p*Delta2p));
    Int[9]  += Fac*((S2xSSP2*SS2)/(64.*Delta1*Delta2p*Delta2p)+(S2xSSP2xAtxS2*SS2)/(64.*Den)-(SS2*SSP2)/(64.*Delta1*Delta2p*Delta2p)-(SS2*SSP2xAtxS2)/(64.*Den)-(S2xSSxSSP*SSxSSP)/(32.*Delta1*Delta2p*Delta2p)-(S2xSSxSSPxAtxS2*SSxSSP)/(32.*Den)+SSxSSP*SSxSSP/(64.*Delta2p*Delta2p)-(AtxS2*SSxSSP*SSxSSP)/(64.*Den)-SSxSSP*SSxSSP/(32.*Delta1*Delta2p*Delta2p)+(S2xAtxS2*SSxSSP*SSxSSP)/(64.*Den)+(SSxSSP*SSxSSPxAtxS2)/(32.*Den)-(SSxSSP*SSxSSPxAtxSS2)/(128.*Den)+(SS2*SSxSSPxAtxSSxSSP)/(256.*Den));
    Int[10] += Fac*(-(S2xSSP2*SS2)/(16.*Delta1*Delta2p*Delta2p)-(S2xSSP2xAtxS2*SS2)/(16.*Den)-(S2xSS2*SSP2)/(16.*Delta1*Delta2p*Delta2p)-(S2xSS2xAtxS2*SSP2)/(16.*Den)+(SS2*SSP2)/(16.*Delta2p*Delta2p)-(AtxS2*SS2*SSP2)/(16.*Den)-(SS2*SSP2)/(8.*Delta1*Delta2p*Delta2p)+(S2xAtxS2*SS2*SSP2)/(16.*Den)+(SS2xAtxS2*SSP2)/(16.*Den)-(SS2xAtxSS2*SSP2)/(64.*Den)+(SS2*SSP2xAtxS2)/(16.*Den)+(S2xSSxSSP*SSxSSP)/(8.*Delta1*Delta2p*Delta2p)+(S2xSSxSSPxAtxS2*SSxSSP)/(8.*Den)-SSxSSP*SSxSSP/(16.*Delta2p*Delta2p)+(AtxS2*SSxSSP*SSxSSP)/(16.*Den)+SSxSSP*SSxSSP/(8.*Delta1*Delta2p*Delta2p)-(S2xAtxS2*SSxSSP*SSxSSP)/(16.*Den)-(SSxSSP*SSxSSPxAtxS2)/(8.*Den)+(SSxSSP*SSxSSPxAtxSS2)/(32.*Den)-(SS2*SSxSSPxAtxSSxSSP)/(64.*Den));
    Int[11] += Fac*((Delta1p*Delta1p*S2xSS2)/(4.*Delta1*Delta2p*Delta2p)+(Delta1p*Delta1p*S2xSS2xAtxS2)/(4.*Den)-(Delta1p*Delta1p*SS2)/(4.*Delta2p*Delta2p)+(AtxS2*Delta1p*Delta1p*SS2)/(4.*Den)+(3*Delta1p*Delta1p*SS2)/(4.*Delta1*Delta2p*Delta2p)-(Delta1p*Delta1p*S2xAtxS2*SS2)/(4.*Den)-(Delta1p*Delta1p*SS2xAtxS2)/(4.*Den)+(Delta1p*Delta1p*SS2xAtxSS2)/(16.*Den));

    /* Vertex Integrals a->0 */
    Inta0[0]  += Fac*(-((16.0*Delta1*Delta1 - Delta2*Delta2)*SS2*(-4.0*Delta1 + SS2))/(256.0*Delta1*Delta1*Delta1*Delta1*Delta2*Delta2));
    Inta0[1]  += Fac*(((-1.0/(16.*Delta1*Delta1) + 1.0/(Delta2*Delta2))*SS2*SS2)/(16.0*Delta1*Delta1));
    Inta0[2]  += Fac*(-(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.0*Delta1*Delta1*Delta2*Delta2));
    Inta0[3]  += Fac*((2.0*S2*S2 + SS2)/(2.0*Delta2*Delta2));
    Inta0[4]  += Fac*((-4.0*AtxS2 - 16.0*S2 + 4.0*S2*S2 + 4.0*S2xAtxS2 + SS2)/(16.0*Delta2*Delta2));
    // Inta0[5]  += 0.0;
    Inta0[6]  += Fac*((4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 4.0*Delta1*SS2 - 4.0*Delta1*Delta1*SS2 - 4.0*S2xAtxS2*SS2 - SS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.0*Delta1*Delta1*Delta2*Delta2));
    Inta0[7]  += Fac*((SS2*(-4.0*Delta1 + SS2))/(16.0*Delta1*Delta1*Delta2*Delta2));
    Inta0[8]  += Fac*((4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 12.0*Delta1*SS2 - 4.0*Delta1*S2*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(32.0*Delta1*Delta2*Delta2));
    Inta0[9]  += Fac*(-(SS2*(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 12*Delta1*SS2 - 4.0*Delta1*Delta1*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2))/(256.0*Delta1*Delta1*Delta2*Delta2));
    // Inta0[10] += 0.0;
    Inta0[11] += Fac*((4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 - 4.0*S2*S2*SS2 - 4.0*S2xAtxS2*SS2 + 4.0*S2*(S2xSS2 + 3*SS2) - 4.0*SS2xAtxS2 + SS2xAtxSS2)/(16.*Delta2*Delta2));

    IntS[0]  += Fac*(-(2.0*CPxSxSSP*PxCPxS + Delta1*PxSSP)/(4.0*Delta1*Delta1*Delta2p*P2));
    IntS[1]  += Fac*((2.0*CPxSxSSP*PxCPxS - Delta1*PxSSP)/(4.0*Delta1*Delta1*Delta2p*P2));
    IntS[2]  += Fac*((32.0*Delta1*Delta1p*PxCPxSP + 32*Delta1p*PxCPxSPxAtxS2 - 16.0*PxCPxSxAtxCPxSxSSP - 32*Delta1p*PxCPxSxAtxSxSP + AtxS2*Delta2p*PxSINP + 3*Delta1*Delta2p*PxSINP - 16.0*Delta1*PxSP2xSSP - 16.0*PxSP2xSSPxAtxS2 - 8.0*AtxS2*PxSSP + 16.0*PxSSPxAtxS2 - 8.0*PxSSP*S2xSP2 - 32*Delta1p*PxCPxS*SxSP + 8.0*PxSSP*SxSP*SxSP + 8.0*PxSSP*SxSPxAtxSxSP)/(32.*Delta1*Delta1*Delta2p*P2));
    IntS[3]  += Fac*((Delta1*Delta2p*PxSINP + 8.0*PxSSP*S2xSP2 + 32*Delta1p*PxCPxS*SxSP - 8.0*PxSSP*SxSP*SxSP)/(32.*Delta1*Delta1*Delta2p*P2));
    IntS[4]  += Fac*((4.0*Delta1*Delta1p*PxCxCPxSS + 4.0*Delta1p*PxCxCPxSSxAtxS2 + 2*CxSPxSSPxAtxS2*PxSS + 2*CxSPxSSP*Delta1*PxSS - 4.0*CxCPxAtxS2*Delta1p*PxSS - 4.0*CxCP*Delta1*Delta1p*PxSS - 2*Delta1p*PxSSxAtxCPxSxSS - 2*Delta1*PxCxSP*SSxSSP - 2*PxCxSPxAtxS2*SSxSSP + PxSSxAtxSxSP*SSxSSP - PxSS*SSxSSPxAtxSxSP)/(8.*Delta1*Delta1*Delta2p*P2));
    IntS[5]  += Fac*((-8.0*Delta1*PxS2xSSP*SS2 - 8.0*PxS2xSSPxAtxS2*SS2 + 8.0*PxSSPxAtxS2*SS2 - 2.0*PxSSxAtxSSxSSP*SS2 - PxSSP*(4.0*Delta1*S2xSS2 + 4.0*S2xSS2xAtxS2 + 4.0*AtxS2*SS2 + 4.0*Delta1*SS2 - 4.0*Delta1*Delta1*SS2 - 4.0*S2xAtxS2*SS2 - 4.0*SS2xAtxS2 + SS2xAtxSS2) + 8.0*Delta1*PxS2xSS*SSxSSP + 8.0*PxS2xSSxAtxS2*SSxSSP - 8.0*PxSSxAtxS2*SSxSSP + 2*PxSSxAtxSS2*SSxSSP + 2*PxSS*(4.0*Delta1*S2xSSxSSP + 4.0*S2xSSxSSPxAtxS2 + 4.0*AtxS2*SSxSSP + 8.0*Delta1*SSxSSP - 4.0*Delta1*Delta1*SSxSSP - 4.0*S2xAtxS2*SSxSSP - 4.0*SSxSSPxAtxS2 + SSxSSPxAtxSS2))/(64.*Delta1*Delta1*Delta2p*P2));
    // IntS[6]  += 0.0;
    // IntS[7]  += 0.0;
    // IntS[8]  += 0.0;
    // IntS[9]  += 0.0;

    IntSa0[0]  += Fac*((Delta1*S2xSS2*((0.125*Delta1 - 0.03125*Delta2)*Delta2 - 0.125*Delta1*SS2) + SS2*(Delta1*Delta2*(-0.0625*Delta1 + 0.0625*Delta1*Delta1 + 0.015625*Delta2 - 0.015625*Delta1*Delta2) + (0.0625*Delta1*Delta1 + 0.125*Delta1*Delta1*Delta1 - 0.00390625*Delta2*Delta2)*SS2))/(Delta1*Delta1*Delta1*Delta1*Delta2*Delta2));
    IntSa0[1]  += Fac*((Delta1*S2xSS2*((-0.125*Delta1 + 0.03125*Delta2)*Delta2 + 0.125*Delta1*SS2) + SS2*(Delta1*Delta1*(-0.0625*Delta1 + 0.015625*Delta2)*Delta2 + (-0.0625*Delta1*Delta1 - 0.125*Delta1*Delta1*Delta1 + 0.00390625*Delta2*Delta2)*SS2))/(Delta1*Delta1*Delta1*Delta1*Delta2*Delta2));
    IntSa0[2]  += Fac*((8.0*CCxS2xSS2xAtxS2 + 8.0*CCxS2xSS2*Delta1 - 8.0*Delta1*Delta2 + 48.0*Delta1*Delta1*Delta2 + 3*Delta1*Delta2*Delta2 - 12*Delta1*Delta1*Delta2*S2 + 4.0*Delta2*S2xAtxS2 - 12*Delta1*Delta2*S2xAtxS2 - 8.0*Delta1*S2xSS2 + 24.0*Delta1*Delta1*S2xSS2 - 2*Delta2*S2xSS2 + 8.0*S2xAtxS2*S2xSS2 + 16.0*S2xSS2xAtxS2 + 16.0*Delta1*S2xSS2xAtxS2 - 4.0*S2xSS2xAtxSS2 + 4.0*Delta1*SS2 + 4.0*Delta1*Delta1*SS2 + 2*Delta2*SS2 - 7*Delta1*Delta2*SS2 - 8.0*Delta1*Delta1*S2*SS2 - 4.0*S2xAtxS2*SS2 - 8.0*Delta1*S2xAtxS2*SS2 + 2*S2xSS2*SS2 - SS2*SS2 - 2*Delta1*SS2*SS2 + AtxS2*(Delta2*(-4.0 + 12*Delta1 + Delta2) - 8.0*S2xSS2 + (4.0 + 8.0*Delta1)*SS2)- 8.0*SS2xAtxS2 - 16.0*Delta1*SS2xAtxS2 - 2*Delta2*SS2xAtxS2 + 2*SS2xAtxSS2 + 4.0*Delta1*SS2xAtxSS2)/(32.*Delta1*Delta1*Delta2*Delta2));
    IntSa0[3]  += Fac*((Delta1*Delta2*(8 + Delta2 - 4.0*Delta1*(3 + S2)) + S2xSS2*(8.0*Delta1 + 8.0*Delta1*Delta1 + 2.0*Delta2 - 2.0*SS2) - (Delta1*(4.0 - 3.0*Delta2) + 2.0*Delta2 + 4.0*Delta1*Delta1*(3 + 2*S2))*SS2 + (1 + 2*Delta1)*SS2*SS2)/(32.*Delta1*Delta1*Delta2*Delta2));
    IntSa0[4]  += Fac*((8.0*CCxS2xSS2xAtxS2*Delta1 + 8.0*CCxS2xSS2*Delta1*Delta1-48.0*Delta1*Delta1*S2xSS2-4.0*Delta1*Delta2*S2xSS2+32*Delta1*Delta1*S2*S2xSS2+16.0*Delta1*S2xAtxS2*S2xSS2+16.0*Delta1*S2xSS2xAtxS2+16.0*Delta1*Delta1*S2xSS2xAtxS2-4.0*Delta2*S2xSS2xAtxS2-4.0*Delta1*S2xSS2xAtxSS2+24.0*Delta1*Delta1*SS2-12*Delta1*Delta2*SS2+4.0*Delta1*Delta1*Delta2*SS2+40*Delta1*Delta1*S2*SS2-16.0*Delta1*Delta1*S2*S2*SS2-8.0*Delta1*S2xAtxS2*SS2-16.0*Delta1*Delta1*S2xAtxS2*SS2+4.0*Delta2*S2xAtxS2*SS2-4.0*AtxS2*(4.0*Delta1*S2xSS2+(-2*Delta1-4.0*Delta1*Delta1+Delta2)*SS2)-8.0*Delta1*SS2xAtxS2-16.0*Delta1*Delta1*SS2xAtxS2+4.0*Delta2*SS2xAtxS2+2*Delta1*SS2xAtxSS2+4.0*Delta1*Delta1*SS2xAtxSS2-Delta2*SS2xAtxSS2)/(32.*Delta1*Delta1*Delta2*Delta2));
    IntSa0[5]  += Fac*((8.0*CCxS2xSS2xAtxS2*Delta2 + 8.0*CCxS2xSS2*Delta1*Delta2+8.0*Delta1*S2xSS2*S2xSS2-4.0*Delta2*S2xSS2xAtxSS2+4.0*AtxS2*Delta1*Delta2*SS2+12*Delta1*Delta1*Delta2*SS2-4.0*Delta1*Delta1*Delta2*S2*SS2-4.0*Delta1*Delta2*S2xAtxS2*SS2-4.0*AtxS2*SS2*SS2-12*Delta1*SS2*SS2-8.0*AtxS2*Delta1*SS2*SS2-20*Delta1*Delta1*SS2*SS2-2*Delta1*Delta2*SS2*SS2+8.0*Delta1*Delta1*S2*SS2*SS2+4.0*S2xAtxS2*SS2*SS2+8.0*Delta1*S2xAtxS2*SS2*SS2-4.0*S2xSS2xAtxS2*(-((2.0+Delta1)*Delta2)+SS2+2.0*Delta1*SS2)-4.0*Delta1*Delta2*SS2xAtxS2+4.0*SS2*SS2xAtxS2+8.0*Delta1*SS2*SS2xAtxS2-2*Delta2*SS2*SS2xAtxS2+Delta1*Delta2*SS2xAtxSS2-SS2*SS2xAtxSS2-2*Delta1*SS2*SS2xAtxSS2+2*S2xSS2*(-4.0*AtxS2*Delta2-12*Delta1*Delta2+6.0*Delta1*Delta1*Delta2+4.0*Delta2*S2xAtxS2+4.0*S2xSS2xAtxS2+4.0*AtxS2*SS2+10*Delta1*SS2-8.0*Delta1*Delta1*SS2-4.0*S2xAtxS2*SS2-4.0*SS2xAtxS2+SS2xAtxSS2))/(128.*Delta1*Delta1*Delta2*Delta2));
    // IntSa0[6]  += 0.0;
    // IntSa0[7]  += 0.0;
    // IntSa0[8]  += 0.0;
    // IntSa0[9]  += 0.0;

    if(i1==0 and i2==0 and i3==0 and i4==1){
      cout<<CCxS2xSS2xAtxS2<<endl;
      exit(0);
    }



  }

  // for(int i1=0; i1<L; i1++)
  // for(int i2=0; i2<L; i2++)
  // for(int i3=0; i3<L; i3++)
  // for(int i4=0; i4<T; i4++)
  // {
  //   int i[4] = {i1,i2,i3,i4};
  //
  //
  // }

  // double P2 = ap.squaredNorm();
  // double L2 = log(P2);
  //
  //


  // vector<double> k(4), kn(4), kp(4), J(4);
  // vector<double> Ss(4)
  //
  //
  // for(int i1=0; i1<L; i1++)
  // for(int i2=0; i2<L; i2++)
  // for(int i3=0; i3<L; i3++)
  // for(int i4=0; i4<T; i4++)
  // {
  //   int i[] = {i1,i2,i3,i4};
  //
  //   for(int mu=0; mu<4; mu++)
  //   {
  //     k[mu]  = 2.0*M_PI/V[mu]*(i[mu]+0.5) - al*sin(2.0*M_PI/V[mu]*(i[mu]+0.5));
  //     kn[mu] = k[mu] + 2.0*ap[mu];
  //     kp[mu] = k[mu] + ap[mu];
  //
  //     J[mu]  = 1.0 - al*cos(2.0*M_PI/V[mu]*(i[mu]+0.5));
  //
  //     S[mu]    = sin(k[mu]/2.0);
  //     SS[mu]   = sin(k[mu]);
  //
  //     S2[mu]   = pow(sin(k[mu]/2.0),2.0);
  //     SS2[mu]  = pow(sin(k[mu]),2.0);
  //
  //     C[mu]    = cos(k[mu]/2.0);
  //     CC[mu]   = cos(k[mu]/2.0);
  //
  //     SP[mu]   = sin(kn[mu]/2.0);
  //     CP[mu]   = cos(kn[mu]/2.0);
  //
  //     SP2[mu]  = pow(sin(kn[mu]/2.0),2.0);
  //     SQ2[mu]  = pow(sin(kp[mu]/2.0),2.0);
  //
  //     SSP[mu]  = sin(kp[mu]);
  //     SSP2[mu] = pow(sin(kp[mu]),2.0);
  //   }
  //
  //   double Jac=0.0, k2=0.0, k4=0.0, k6=0.0, k10p = 0.0;
  //
  //   for(int mu=0; mu<4; mu++)
  //   {
  //     Jac *= J[mu];
  //
  //     k2   += 4.0*S2[mu];
  //     k4   += pow(4.0*S2[mu],2.0);
  //     k6   += pow(4.0*S2[mu],3.0);
  //     k10p += 1244.0*pow(S2[mu],2.0)*S2[(mu+1)%3]*S2[(mu+2)%3]*S2[(mu+3)%3];
  //   }
  //
  //   for(int mu=0; mu<4; mu++){
  //     for(int mu=0; mu<4; mu++){
  //       A[mu][nu] = (c1*k2*k2*(k2 - k2t12) - 0.5*c1*c1*k2*(k2*k2*k2 - 2.0*k2t12*k4 + k2*(k4 - 2.0*k4p12) + 2.0*k6) + c1*c1*c1*(0.5*k4.0*(k2*k2*k2 - k2*k4 + 2.0*k6) + 4.0*k10p))/Del;
  //     }
  //   }
  //
  //
  //
  //   double Del = (k2 - c1*k4)*(k2 - c1*(k2*k2 + k4) + c1*c1*0.5*(k2*k2*k2 + 2.0*k6 - k2*k4)) - 4.0*c1*c1*c1*k10p;
  //
  //
  //
  //
  // }







  // }


  //
  // valarray<v1d_t> k(2),J(2),Ss(2),Ss2(2),SSs(2),SSs2(2),C(2),CC(2);
  // valarray<v1d_t> Sp(4),Cp(4),Sp2(4),Sq2(4),SSp(4),SSp2(4);
  //
  // for(int i=0; i<2; i++)
  // {
  //   k[i]    = vd_t(0.0,dim[i]);
  //   J[i]    = vd_t(0.0,dim[i]);
  //   Ss[i]   = vd_t(0.0,dim[i]);
  //   Ss2[i]  = vd_t(0.0,dim[i]);
  //   SSs[i]  = vd_t(0.0,dim[i]);
  //   SSs2[i] = vd_t(0.0,dim[i]);
  //   C[i]    = vd_t(0.0,dim[i]);
  //   CC[i]   = vd_t(0.0,dim[i]);
  // }
  //
  // #pragma omp parallel for
  // for(int i=0; i<2; i++)
  // for(int j=0; j<dim[i]; j++)
  // {
  //   k[i][j] = 2.0*M_PI/dim[i]*(j+0.5) - al*sin(2.0*M_PI/dim[i]*(j+0.5));
  //   J[i][j] = 1.0 - al*cos(2.0*M_PI/dim[i]*(j+0.5));
  // }
  //
  // #pragma omp parallel for
  // for(int i=0; i<2; i++)
  // {
  //   Ss[i]   = sin(k[i]/2.0);
  //   Ss2[i]  = Ss[i]*Ss[i];
  //   SSs[i]  = sin(k[i]);
  //   SSs2[i] = SSs[i]*SSs[i];
  //   C[i]    = cos(k[i]/2.0);
  //   CC[i]   = cos(k[i]);
  // }
  //
  // #pragma omp parallel for
  // for(int i=0; i<4; i++)
  // {
  //   // int j = (int)(i/3);
  //
  //   Sp[i]   = sin((k[(int)(i/3)] + 2.0*ap[i])/2.0);
  //   Cp[i]   = cos((k[(int)(i/3)] + 2.0*ap[i])/2.0);
  //   Sp2[i]  = Sp[i]*Sp[i];
  //   Sq2[i]  = sin((k[(int)(i/3)] + ap[i])/2.0)*sin((k[(int)(i/3)] + ap[i])/2.0);
  //   SSp[i]  = sin(k[(int)(i/3)] + ap[i]);
  //   SSp2[i] = SSp[i]*SSp[i];
  // }
  //
  // v4d_t Jac(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  //
  // valarray<v4d_t> s   (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> sp  (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> s2  (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> sp2 (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> sq2 (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> ss  (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> ssp (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> ss2 (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> ssp2(v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> co  (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> cp  (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  // valarray<v4d_t> cc  (v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4);
  //
  //
  // #pragma omp parallel for collapse(4)
  // for(int i1=0; i1<L; i1++)
  // for(int i2=0; i2<L; i2++)
  // for(int i3=0; i3<L; i3++)
  // for(int i4=0; i4<T; i4++)
  // {
  //   int ind[]={i1,i2,i3,i4};
  //
  //   Jac[i1][i2][i3][i4] = J[0][i1]*J[0][i2]*J[0][i3]*J[1][i4];
  //
  //   for(int i=0; i<4; i++)
  //   {
  //     int j=(int)(i/3);
  //
  //     s[i][i1][i2][i3][i4] = Ss[j][ind[i]];
  //     // s[1][i1][i2][i3][i4] = Ss[0][i2];
  //     // s[2][i1][i2][i3][i4] = Ss[0][i3];
  //     // s[3][i1][i2][i3][i4] = Ss[1][i4];
  //     /**/
  //     sp[i][i1][i2][i3][i4] = Sp[i][ind[i]];
  //     // sp[1][i1][i2][i3][i4] = Sp[1][i2];
  //     // sp[2][i1][i2][i3][i4] = Sp[2][i3];
  //     // sp[3][i1][i2][i3][i4] = Sp[3][i4];
  //     /**/
  //     s2[i][i1][i2][i3][i4] = Ss2[j][ind[i]];
  //     // s2[1][i1][i2][i3][i4] = Ss2[0][i2];
  //     // s2[2][i1][i2][i3][i4] = Ss2[0][i3];
  //     // s2[3][i1][i2][i3][i4] = Ss2[1][i4];
  //     /**/
  //     sp2[i][i1][i2][i3][i4] = Sp2[i][ind[i]];
  //     // sp2[1][i1][i2][i3][i4] = Sp2[1][i2];
  //     // sp2[2][i1][i2][i3][i4] = Sp2[2][i3];
  //     // sp2[3][i1][i2][i3][i4] = Sp2[3][i4];
  //     /**/
  //     sq2[i][i1][i2][i3][i4] = Sq2[i][ind[i]];
  //     // sq2[1][i1][i2][i3][i4] = Sq2[1][i2];
  //     // sq2[2][i1][i2][i3][i4] = Sq2[2][i3];
  //     // sq2[3][i1][i2][i3][i4] = Sq2[3][i4];
  //     /**/
  //     ss[i][i1][i2][i3][i4] = SSs[j][ind[i]];
  //     // ss[1][i1][i2][i3][i4] = SSs[0][i2];
  //     // ss[2][i1][i2][i3][i4] = SSs[0][i3];
  //     // ss[3][i1][i2][i3][i4] = SSs[1][i4];
  //     /**/
  //     ssp[i][i1][i2][i3][i4] = SSp[i][ind[i]];
  //     // ssp[1][i1][i2][i3][i4] = SSp[1][i2];
  //     // ssp[2][i1][i2][i3][i4] = SSp[2][i3];
  //     // ssp[3][i1][i2][i3][i4] = SSp[3][i4];
  //     /**/
  //     ss2[i][i1][i2][i3][i4] = SSs2[j][ind[i]];
  //     // ss2[1][i1][i2][i3][i4] = SSs2[0][i2];
  //     // ss2[2][i1][i2][i3][i4] = SSs2[0][i3];
  //     // ss2[3][i1][i2][i3][i4] = SSs2[1][i4];
  //     /**/
  //     ssp2[i][i1][i2][i3][i4] = SSp2[i][ind[i]];
  //     // ssp2[1][i1][i2][i3][i4] = SSp2[1][i2];
  //     // ssp2[2][i1][i2][i3][i4] = SSp2[2][i3];
  //     // ssp2[3][i1][i2][i3][i4] = SSp2[3][i4];
  //     /**/
  //     co[i][i1][i2][i3][i4] = C[j][ind[i]];
  //     // co[1][i1][i2][i3][i4] = C[0][i2];
  //     // co[2][i1][i2][i3][i4] = C[0][i3];
  //     // co[3][i1][i2][i3][i4] = C[1][i4];
  //     /**/
  //     cc[i][i1][i2][i3][i4] = CC[j][ind[i]];
  //     // cc[1][i1][i2][i3][i4] = CC[0][i2];
  //     // cc[2][i1][i2][i3][i4] = CC[0][i3];
  //     // cc[3][i1][i2][i3][i4] = CC[1][i4];
  //     /**/
  //     cp[i][i1][i2][i3][i4] = Cp[i][ind[i]];
  //     // cp[1][i1][i2][i3][i4] = Cp[1][i2];
  //     // cp[2][i1][i2][i3][i4] = Cp[2][i3];
  //     // cp[3][i1][i2][i3][i4] = Cp[3][i4];
  //   }
  // }
  //
  // /* Gluon Propagator Symanzik-improved */
  //
  // v4d_t k2 = 4.0*(s2[0]+s2[1]+s2[2]+s2[3]);
  // v4d_t k4 = 16.0*(s2[0]*s2[0]+s2[1]*s2[1]+s2[2]*s2[2]+s2[3]*s2[3]);
  // v4d_t k6 = 64.0*(s2[0]*s2[0]*s2[0]+s2[1]*s2[1]*s2[1]+s2[2]*s2[2]*s2[2]+s2[3]*s2[3]*s2[3]);
  // v4d_t k10p = 1244.0*(s2[0]*s2[0]*s2[1]*s2[2]*s2[3]+s2[0]*s2[1]*s2[1]*s2[2]*s2[3]+s2[0]*s2[1]*s2[2]*s2[2]*s2[3]+s2[0]*s2[1]*s2[2]*s2[3]*s2[3]);
  //
  // v4d_t k2t12 = 4.0*(s2[2]+s2[3]);
  // v4d_t k2t14 = 4.0*(s2[1]+s2[2]);
  // v4d_t k4p12 = 16.0*(s2[2]*s2[3]);
  // v4d_t k4p14 = 16.0*(s2[1]*s2[2]);
  //
  // v4d_t Del = (k2 - c1*k4)*(k2 - c1*(k2*k2 + k4) + c1*c1*0.5*(k2*k2*k2 + 2.0*k6 - k2*k4)) - 4.0*c1*c1*c1*k10p;
  //
  //
  // valarray<valarray<v4d_t>> A(valarray<v4d_t>(v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),4),4);
  //
  // A[0][1] = (c1*k2*k2*(k2 - k2t12) - 0.5*c1*c1*k2*(k2*k2*k2 - 2.0*k2t12*k4 + k2*(k4 - 2.0*k4p12) + 2.0*k6) + c1*c1*c1*(0.5*k4.0*(k2*k2*k2 - k2*k4 + 2.0*k6) + 4.0*k10p))/Del;
  // A[0][3] = (c1*k2*k2*(k2 - k2t14) - 0.5*c1*c1*k2*(k2*k2*k2 - 2.0*k2t14.0*k4 + k2*(k4 - 2.0*k4p14) + 2.0*k6) + c1*c1*c1*(0.5*k4.0*(k2*k2*k2 - k2*k4 + 2.0*k6) + 4.0*k10p))/Del;
  //
  // #pragma omp parallel for collapse(4)
  // for(int i1=0; i1<L; i1++)
  // for(int i2=0; i2<L; i2++)
  // for(int i3=0; i3<L; i3++)
  // for(int i4=0; i4<T; i4++)
  // {
  //   A[0][2][i1][i2][i3][i4] = A[0][1][i1][i3][i2][i4];
  //   A[1][2][i1][i2][i3][i4] = A[0][1][i3][i2][i1][i4];
  //   A[1][3][i1][i2][i3][i4] = A[0][3][i2][i1][i3][i4];
  //   A[2][3][i1][i2][i3][i4] = A[0][3][i3][i2][i1][i4];
  // }
  //
  // A[1][0] = A[0][1];
  // A[3][0] = A[0][3];
  // A[2][0] = A[0][2];
  // A[2][1] = A[1][2];
  // A[3][1] = A[1][3];
  // A[3][2] = A[2][3];
  //
  // // for(int i=0;i<4;i++)
  // // {
  // //   for(int j=0;j<4;j++)
  // //   {
  // //     cout<<A[i][j][0][0][0][0]<<" ";
  // //   }
  // //   cout<<endl;
  // // }
  // // cout<<endl;
  // //
  // // for(int i=0;i<4;i++)
  // // {
  // //   for(int j=0;j<4;j++)
  // //   {
  // //     cout<<A[i][j][0][0][0][1]<<" ";
  // //   }
  // //   cout<<endl;
  // // }
  //
  //
  //
  // /* Invariants */
  //
  // /* Without gluon propagator */
  // v4d_t S2  (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SP2 (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SQ2 (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SS2 (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SSP2(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  //
  // v4d_t SxSP       (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CxCP       (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SSxSSP     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xSP2     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xSS2     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xSSP2    (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SP2xSSP2   (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xSSxSSP  (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CPxSxSS    (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CPxSxSSP   (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CPxSPxSSP  (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CxCPxSSxSSP(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CxSPxSSP   (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CCxS2xSS2  (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  //
  // v4d_t PxSS     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSSP    (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxCPxS   (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxS2xSSP (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxS2xSS  (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxCxCPxSS(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxCPxSP  (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxCxSP   (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSP2xSSP(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  //
  // /* With gluon propagator */
  // v4d_t CCxS2xSS2xAtxS2(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xSS2xAtxSS2  (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xSS2xAtxS2   (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SS2xAtxSS2     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SS2xAtxS2      (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xAtxS2       (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t AtxS2          (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  //
  // v4d_t CPxSxSSPxAtxCPxSxSSP(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CPxSxSSxAtxSSxSSP   (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CPxSxSSPxAtxSxSP    (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CPxSPxSSPxAtxS2     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CxCPxSSxSSPxAtxS2   (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CxSPxSSPxAtxS2      (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t CxCPxAtxS2          (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xSSxSSPxAtxS2     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t S2xSSP2xAtxS2       (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SP2xSSP2xAtxS2      (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SSxSSPxAtxSSxSSP    (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SSxSSPxAtxSxSP      (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SSxSSPxAtxSS2       (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SSxSSPxAtxS2        (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SxSPxAtxSxSP        (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SSP2xAtxS2          (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t SP2xAtxS2           (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  //
  // v4d_t PxCPxSxAtxCPxSxSSP (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxCPxSxAtxSxSP     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxCPxSPxAtxS2      (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxCxCPxSSxAtxS2    (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxCxSPxAtxS2       (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSSPxAtxS2        (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSP2xSSPxAtxS2    (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxS2xSSPxAtxS2     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxS2xSSxAtxS2      (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSSxAtxCPxSxSS    (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSSxAtxSSxSSP     (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSSxAtxSxSP       (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSSxAtxSS2        (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  // v4d_t PxSSxAtxS2         (v3d_t(v2d_t(v1d_t(0.0,T),L),L),L);
  //
  // double COSP = 0.0;
  // // double PxSINP = P2;
  //
  // #pragma omp parallel for
  // for(int i=0; i<4; i++)
  // {
  //   S2   += s2[i];
  //   SP2  += sp2[i];
  //   SQ2  += sq2[i];
  //   SS2  += ss2[i];
  //   SSP2 += ssp2[i];
  //   /**/
  //   SxSP        += s[i]*sp[i];
  //   CxCP        += co[i]*cp[i];
  //   SSxSSP      += ss[i]*ssp[i];
  //   S2xSP2      += s2[i]*sp2[i];
  //   S2xSS2      += s2[i]*ss2[i];
  //   S2xSSP2     += s2[i]*ssp2[i];
  //   SP2xSSP2    += sp2[i]*ssp2[i];
  //   S2xSSxSSP   += s2[i]*ss[i]*ssp[i];
  //   CPxSxSS     += cp[i]*s[i]*ss[i];
  //   CPxSxSSP    += cp[i]*s[i]*ssp[i];
  //   CPxSPxSSP   += cp[i]*sp[i]*ssp[i];
  //   CxCPxSSxSSP += co[i]*cp[i]*ss[i]*ssp[i];
  //   CxSPxSSP    += co[i]*sp[i]*ssp[i];
  //   CCxS2xSS2   += cc[i]*s2[i]*ss2[i];
  //   /**/
  //   COSP += cos(ap[i]);
  //   /**/
  //   if(ap[i])
  //   {
  //     PxSS      += P2/Np0*(1.0/sin(ap[i]))*ss[i];
  //     PxSSP     += P2/Np0*(1.0/sin(ap[i]))*ssp[i];
  //     PxCPxS    += P2/Np0*(1.0/sin(ap[i]))*cp[i]*s[i];
  //     PxS2xSSP  += P2/Np0*(1.0/sin(ap[i]))*s2[i]*ssp[i];
  //     PxS2xSS   += P2/Np0*(1.0/sin(ap[i]))*s2[i]*ss[i];
  //     PxCxCPxSS += P2/Np0*(1.0/sin(ap[i]))*co[i]*cp[i]*ss[i];
  //     PxCPxSP   += P2/Np0*(1.0/sin(ap[i]))*cp[i]*sp[i];
  //     PxCxSP    += P2/Np0*(1.0/sin(ap[i]))*co[i]*sp[i];
  //     PxSP2xSSP += P2/Np0*(1.0/sin(ap[i]))*sp2[i]*ssp[i];
  //   }
  // }
  //
  // // #pragma omp parallel for collapse(2)
  // for(int i=0; i<4; i++)
  // for(int j=0; j<4; j++)
  // {
  //   CCxS2xSS2xAtxS2 += cc[i]*s2[i]*ss2[i]*(A[i][j]*s2[j]);
  //   S2xSS2xAtxSS2   += s2[i]*ss2[i]*(A[i][j]*ss2[j]);
  //   S2xSS2xAtxS2    += s2[i]*ss2[i]*(A[i][j]*s2[j]);
  //   SS2xAtxSS2      += ss2[i]*(A[i][j]*ss2[j]);
  //   SS2xAtxS2       += ss2[i]*(A[i][j]*s2[j]);
  //   S2xAtxS2        += s2[i]*(A[i][j]*s2[j]);
  //   AtxS2           += A[i][j]*s2[j];
  //   /**/
  //   CPxSxSSPxAtxCPxSxSSP += cp[i]*s[i]*ssp[i]*(A[i][j]*cp[j]*s[j]*ssp[j]);
  //   CPxSxSSxAtxSSxSSP    += cp[i]*s[i]*ss[i]*(A[i][j]*ss[j]*ssp[j]);
  //   CPxSxSSPxAtxSxSP     += cp[i]*s[i]*ssp[i]*(A[i][j]*s[j]*sp[j]);
  //   CPxSPxSSPxAtxS2      += cp[i]*sp[i]*ssp[i]*(A[i][j]*s2[j]);
  //   CxCPxSSxSSPxAtxS2    += co[i]*cp[i]*ss[i]*ssp[i]*(A[i][j]*s2[j]);
  //   CxSPxSSPxAtxS2       += co[i]*sp[i]*ssp[i]*(A[i][j]*s2[j]);
  //   CxCPxAtxS2           += co[i]*cp[i]*(A[i][j]*s2[j]);
  //   S2xSSxSSPxAtxS2      += s2[i]*ss[i]*ssp[i]*(A[i][j]*s2[j]);
  //   S2xSSP2xAtxS2        += s2[i]*ssp2[i]*(A[i][j]*s2[j]);
  //   SP2xSSP2xAtxS2       += sp2[i]*ssp2[i]*(A[i][j]*s2[j]);
  //   SSxSSPxAtxSSxSSP     += ss[i]*ssp[i]*(A[i][j]*ss[j]*ssp[j]);
  //   SSxSSPxAtxSxSP       += ss[i]*ssp[i]*(A[i][j]*s[j]*sp[j]);
  //   SSxSSPxAtxSS2        += ss[i]*ssp[i]*(A[i][j]*ss2[j]);
  //   SSxSSPxAtxS2         += ss[i]*ssp[i]*(A[i][j]*s2[j]);
  //   SxSPxAtxSxSP         += s[i]*sp[i]*(A[i][j]*s[j]*sp[j]);
  //   SSP2xAtxS2           += ssp2[i]*(A[i][j]*s2[j]);
  //   SP2xAtxS2            += sp2[i]*(A[i][j]*s2[j]);
  //   /**/
  //   if(ap[i]!=0.0)
  //   {
  //     PxCPxSxAtxCPxSxSSP += P2/Np0*(1.0/sin(ap[i]))*cp[i]*s[i]*(A[i][j]*cp[j]*s[j]*ssp[j]);
  //     PxCPxSxAtxSxSP     += P2/Np0*(1.0/sin(ap[i]))*cp[i]*s[i]*(A[i][j]*s[j]*sp[j]);
  //     PxCPxSPxAtxS2      += P2/Np0*(1.0/sin(ap[i]))*sp[i]*(A[i][j]*s2[j]);
  //     PxCxCPxSSxAtxS2    += P2/Np0*(1.0/sin(ap[i]))*co[i]*cp[i]*ss[i]*(A[i][j]*s2[j]);
  //     PxCxSPxAtxS2       += P2/Np0*(1.0/sin(ap[i]))*co[i]*sp[i]*(A[i][j]*s2[j]);
  //     PxSSPxAtxS2        += P2/Np0*(1.0/sin(ap[i]))*ssp[i]*(A[i][j]*s2[j]);
  //     PxSP2xSSPxAtxS2    += P2/Np0*(1.0/sin(ap[i]))*sp2[i]*ssp[i]*(A[i][j]*s2[j]);
  //     PxS2xSSPxAtxS2     += P2/Np0*(1.0/sin(ap[i]))*s2[i]*ssp[i]*(A[i][j]*s2[j]);
  //     PxS2xSSxAtxS2      += P2/Np0*(1.0/sin(ap[i]))*s2[i]*ss[i]*(A[i][j]*s2[j]);
  //     PxSSxAtxCPxSxSS    += P2/Np0*(1.0/sin(ap[i]))*ss[i]*(A[i][j]*cp[j]*s[j]*ss[j]);
  //     PxSSxAtxSSxSSP     += P2/Np0*(1.0/sin(ap[i]))*ss[i]*(A[i][j]*ss[j]*ssp[j]);
  //     PxSSxAtxSxSP       += P2/Np0*(1.0/sin(ap[i]))*ss[i]*(A[i][j]*s[j]*sp[j]);
  //     PxSSxAtxSS2        += P2/Np0*(1.0/sin(ap[i]))*ss[i]*(A[i][j]*ss2[j]);
  //     PxSSxAtxS2         += P2/Np0*(1.0/sin(ap[i]))*ss[i]*(A[i][j]*s2[j]);
  //   }
  // }
  //
  //
  // /* Denominators */
  // v4d_t Delta1  = S2;
  // v4d_t Delta1p = SQ2;
  // v4d_t Delta2  = SS2 + 4.0*(S2*S2);
  // v4d_t Delta2p = SSP2 + 4.0*(SQ2*SQ2);
  // v4d_t Den = (Delta1*Delta1*Delta2p*Delta2p);
  //
  // /* Integrals */
  // /* Vertex Integrals */
  // valarray<v4d_t> IntX(v4d_t(v3d_t(v2d_t(v1d_t(0.0,T),L),L),L),12);
  // vd_t Int(0.0,12);
  //
  // IntX[0]  = -CPxSxSSP*CPxSxSSP/(4.*Den) + SSP2/(4.*Delta1*Delta2p*Delta2p);
  // IntX[1]  =  CPxSxSSP*CPxSxSSP/(4.*Den);
  // IntX[2]  = -CPxSxSSPxAtxCPxSxSSP/(4.*Den)+(CPxSPxSSPxAtxS2*Delta1p)/(Den)-(CPxSxSSPxAtxSxSP*Delta1p)/(Den)+(CPxSPxSSP*Delta1p)/(Delta1*Delta2p*Delta2p)+(Delta1p*Delta1p*SP2)/(Delta1*Delta2p*Delta2p)+(Delta1p*Delta1p*SP2xAtxS2)/(Den)-SP2xSSP2/(4.*Delta1*Delta2p*Delta2p)-SP2xSSP2xAtxS2/(4.*Den)+SSP2xAtxS2/(4.*Den)-(CPxSxSSP*Delta1p*SxSP)/(Den)-(Delta1p*Delta1p*SxSP*SxSP)/(Den)-(Delta1p*Delta1p*SxSPxAtxSxSP)/(Den);
  // IntX[3]  =  (CPxSxSSP*Delta1p*SxSP)/(Den)+(Delta1p*Delta1p*SxSP*SxSP)/(Den);
  // IntX[4]  =  (CPxSPxSSPxAtxS2*Delta1p)/(4.*Den)-(CPxSxSSPxAtxSxSP*Delta1p)/(4.*Den)+(CPxSPxSSP*Delta1p)/(4.*Delta1*Delta2p*Delta2p)-(AtxS2*Delta1p*Delta1p)/(4.*Den)-(3*Delta1p*Delta1p)/(4.*Delta1*Delta2p*Delta2p)-(Delta1p*Delta1p*S2xSP2)/(4.*Den)+(Delta1p*Delta1p*SP2)/(4.*Delta1*Delta2p*Delta2p)+(Delta1p*Delta1p*SP2xAtxS2)/(4.*Den)-(SP2*SSP2)/(16.*Delta1*Delta2p*Delta2p)-(SP2xAtxS2*SSP2)/(16.*Den)-(CPxSxSSP*Delta1p*SxSP)/(4.*Den)+(SSP2*SxSP*SxSP)/(16.*Den)+(SSP2*SxSPxAtxSxSP)/(16.*Den);
  // IntX[5]  =  -Delta1p*Delta1p/(4.*Delta1*Delta2p*Delta2p) + (Delta1p*Delta1p*S2xSP2)/(4.*Den) + (CPxSxSSP*Delta1p*SxSP)/(4.*Den) - (SSP2*SxSP*SxSP)/(16.*Den);
  // IntX[6]  = CPxSxSSPxAtxCPxSxSSP/(4.*Den)+SP2xSSP2/(4.*Delta1*Delta2p*Delta2p)+SP2xSSP2xAtxS2/(4.*Den)+(AtxS2*SSP2)/(4.*Den)+(S2xSP2*SSP2)/(4.*Den)-(SP2*SSP2)/(4.*Delta1*Delta2p*Delta2p)-(SP2xAtxS2*SSP2)/(4.*Den)-SSP2xAtxS2/(4.*Den);
  // IntX[7]  =  -(S2xSP2*SSP2)/(4.*Den);
  // IntX[8]  = (CPxSxSSxAtxSSxSSP*Delta1p)/(16.*Den)-(CxCPxSSxSSPxAtxS2*Delta1p)/(8.*Den)-(CxCPxSSxSSP*Delta1p)/(8.*Delta1*Delta2p*Delta2p)+(CxCPxAtxS2*Delta1p*SSxSSP)/(8.*Den)+(CxCP*Delta1p*SSxSSP)/(8.*Delta1*Delta2p*Delta2p);
  // IntX[9]  = (S2xSSP2*SS2)/(64.*Delta1*Delta2p*Delta2p)+(S2xSSP2xAtxS2*SS2)/(64.*Den)-(SS2*SSP2)/(64.*Delta1*Delta2p*Delta2p)-(SS2*SSP2xAtxS2)/(64.*Den)-(S2xSSxSSP*SSxSSP)/(32.*Delta1*Delta2p*Delta2p)-(S2xSSxSSPxAtxS2*SSxSSP)/(32.*Den)+SSxSSP*SSxSSP/(64.*Delta2p*Delta2p)-(AtxS2*SSxSSP*SSxSSP)/(64.*Den)-SSxSSP*SSxSSP/(32.*Delta1*Delta2p*Delta2p)+(S2xAtxS2*SSxSSP*SSxSSP)/(64.*Den)+(SSxSSP*SSxSSPxAtxS2)/(32.*Den)-(SSxSSP*SSxSSPxAtxSS2)/(128.*Den)+(SS2*SSxSSPxAtxSSxSSP)/(256.*Den);
  // IntX[10] = -(S2xSSP2*SS2)/(16.*Delta1*Delta2p*Delta2p)-(S2xSSP2xAtxS2*SS2)/(16.*Den)-(S2xSS2*SSP2)/(16.*Delta1*Delta2p*Delta2p)-(S2xSS2xAtxS2*SSP2)/(16.*Den)+(SS2*SSP2)/(16.*Delta2p*Delta2p)-(AtxS2*SS2*SSP2)/(16.*Den)-(SS2*SSP2)/(8.*Delta1*Delta2p*Delta2p)+(S2xAtxS2*SS2*SSP2)/(16.*Den)+(SS2xAtxS2*SSP2)/(16.*Den)-(SS2xAtxSS2*SSP2)/(64.*Den)+(SS2*SSP2xAtxS2)/(16.*Den)+(S2xSSxSSP*SSxSSP)/(8.*Delta1*Delta2p*Delta2p)+(S2xSSxSSPxAtxS2*SSxSSP)/(8.*Den)-SSxSSP*SSxSSP/(16.*Delta2p*Delta2p)+(AtxS2*SSxSSP*SSxSSP)/(16.*Den)+SSxSSP*SSxSSP/(8.*Delta1*Delta2p*Delta2p)-(S2xAtxS2*SSxSSP*SSxSSP)/(16.*Den)-(SSxSSP*SSxSSPxAtxS2)/(8.*Den)+(SSxSSP*SSxSSPxAtxSS2)/(32.*Den)-(SS2*SSxSSPxAtxSSxSSP)/(64.*Den);
  // IntX[11] =  (Delta1p*Delta1p*S2xSS2)/(4.*Delta1*Delta2p*Delta2p)+(Delta1p*Delta1p*S2xSS2xAtxS2)/(4.*Den)-(Delta1p*Delta1p*SS2)/(4.*Delta2p*Delta2p)+(AtxS2*Delta1p*Delta1p*SS2)/(4.*Den)+(3*Delta1p*Delta1p*SS2)/(4.*Delta1*Delta2p*Delta2p)-(Delta1p*Delta1p*S2xAtxS2*SS2)/(4.*Den)-(Delta1p*Delta1p*SS2xAtxS2)/(4.*Den)+(Delta1p*Delta1p*SS2xAtxSS2)/(16.*Den);
  //
  // #pragma omp parallel for
  // for(int i=0; i<12; i++)
  // for(int i1=0; i1<L; i1++)
  // for(int i2=0; i2<L; i2++)
  // for(int i3=0; i3<L; i3++)
  // for(int i4=0; i4<T; i4++)
  // Int[i] += (IntX[i][i1][i2][i3][i4])*(Jac[i1][i2][i3][i4])*16.0*M_PI*M_PI/(L*L*L*T);
  //
  // for(int i=0;i<12;i++) cout<<Int[i]<<endl;


  exit(0);
}
