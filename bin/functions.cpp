
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include "global.hpp"
// #include "aliases.hpp"

// using namespace Eigen;
using namespace std;


vector<double> compute_V(double *Int)
{
  vector<double> V(5);
  double VV00, VV10, VV11, VV12;
  double csw2 = csw*csw;

  double h1[5] = {  4.0, -4.0, -2.0,  2.0, 0.0};
  double h2[5] = { 16.0, 16.0,  4.0,  4.0, 0.0};
  double h3[5] = {  1.0, -1.0,  0.5, -0.5, 0.0};
  double h4[5];
  for (int i=0; i<5; i++){
    h4[i] = (h2[i]/4.0-1.0)/3.0;

    VV00 = (1.0+(2.0+alpha)*h4[i])*Int[0] + (alpha+3*h4[i])*Int[1];
    VV10 = Int[2]+alpha*Int[3]+(Int[4]+alpha*Int[5])*h1[i] + (Int[6]+alpha*Int[7])*h4[i];
    VV11 = Int[8]*(h1[i]+4.0*h4[i]);
    VV12 = Int[9]*h1[i] + Int[10]*h3[i] + Int[11]*h4[i];

    V[i] = VV00 + VV10 + csw*VV11 + csw2*VV12;
  }


  return V;
}

double compute_S(double *Int)
{
  double S = Int[0] + alpha*Int[1] + Int[2] + alpha*Int[3] + csw*Int[4] + csw*csw*Int[5];

  return S ;
}

vector<double> compute_Z(double *Int, double *IntS)
{
  vector<double> Z(6);

  vector<double> V=compute_V(Int);
  double S=compute_S(IntS);

  double g02 = 6.0/beta;
  double as = g02/(16.0*M_PI*M_PI)*4.0/3.0;

  Z[0] = -S*as;
  for(int i=0; i<5; i++)
  {
    Z[i+1] = Z[0] - V[i]*as;
  }
  return Z;
}

inline double norm3(const array<double,4>& mom)
{
  return mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2];
}

void find_eqmoms(){

  array<double,4> shift={0.0,0.0,0.0,0.0};
  if(APBC) shift[0]+=0.5;

  ifstream input("mom_list.txt");
  if(!input.good())
  {
    cerr<<"Error opening \"mom_list.txt\"."<<endl;
    exit(1);
  }

  while(!input.eof())
  {
    array<double,4> np,ap;
    /* np = {T,L0,L1,L2} */
    /* ap = {L0,L1,L2,T} */

    for(int mu=0;mu<4;mu++){
      input>>np[mu];
      ap[(mu+3)%4] = 2.0*M_PI*(np[mu]+shift[mu])/V[(mu+3)%4];
    }
    if(input.good())
    {
      ap_list.push_back(ap);
    }
  }

  // Find equivalent moms
  ap_eq.push_back(ap_list[0]);
  tag_list.push_back(0);

  moms=ap_list.size();
  int tag=0, tag_aux=0;
  double eps=1.0e-15;

  for(int imom=0; imom<moms; imom++)
  {
    int count_no = 0;

    for(int j=0;j<imom;j++)
    {
      /* H(3,1) condition */
      bool cond{
        fabs(fabs(ap_list[j][3])-fabs(ap_list[imom][3]))<eps*(fabs(ap_list[j][3])+fabs(ap_list[imom][3])) &&
        fabs(norm3(ap_list[j])-norm3(ap_list[imom]))<eps*(norm3(ap_list[j])+norm3(ap_list[imom]))
      };

      /* identity condition */
      // bool cond{ap_list[j][0]==ap_list[imom][0] && ap_list[j][1]==ap_list[imom][1] && ap_list[j][2]==ap_list[imom][2] && ap_list[j][3]==ap_list[imom][3]};  /*test*/

      if(cond)  tag_aux = tag_list[j];
      else      count_no++;

      if(count_no==imom)
      {
          tag++;
          tag_list.push_back(tag);
          ap_eq.push_back(ap_list[imom]);
      }
      else if(j==imom-1)
      {
        tag_list.push_back(tag_aux);
      }
    }
  }
  eqmoms = ap_eq.size();
}
