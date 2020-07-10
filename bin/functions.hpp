#ifndef FUNC_HPP
#define FUNC_HPP


using namespace Eigen;
using namespace std;

// inline Vector4d cwiseProduct(Vector4d a, Vector4d b)
// {
//   return a.cwiseProduct(b);
// }
// inline Vector4d cwiseProduct(Vector4d a, Vector4d b, Vector4d c)
// {
//   return a.cwiseProduct(b).cwiseProduct(c);
// }
// inline Vector4d cwiseProduct(Vector4d a, Vector4d b, Vector4d c, Vector4d d)
// {
//   return a.cwiseProduct(b).cwiseProduct(c).cwiseProduct(d);
// }
// inline Vector4d cwiseProduct(Vector4d a, Vector4d b, Vector4d c, Vector4d d, Vector4d e)
// {
//   return a.cwiseProduct(b).cwiseProduct(c).cwiseProduct(d).cwiseProduct(e);
// }

// inline double prod(Vector4d a)
// {
//   return a.dot(Id);
// }
// inline double prod(Vector4d a,Vector4d b)
// {
//   return a.dot(b);
// };
// inline double prod(Vector4d a,Vector4d b,Vector4d c)
// {
//   return cwiseProduct(a,b,c).dot(Id);
// };
// inline double prod(Vector4d a,Vector4d b,Vector4d c,Vector4d d)
// {
//   return cwiseProduct(a,b,c,d).dot(Id);
// };
// inline double prod(Vector4d a,Vector4d b,Vector4d c,Vector4d d,Vector4d e)
// {
//   return cwiseProduct(a,b,c,d,e).dot(Id);
// };

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

      ktpo2[mu][nu]=1.0;
      ktso2[mu][nu]=0.0;

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

  T Deltakt=(kt2-c1*kt4)*(kt2-c1*(kt22+kt4)+0.5*c12*(kt23+2.0*kt6-kt2*kt4));
  for(int rho=0;rho<4;rho++) Deltakt-=4.0*c13*kt4_dir[rho]*ktpo2[rho][rho];

  for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++){
      A[mu][nu]=(1.0 - kron_delta[mu][nu])/Deltakt*(kt22-c1*kt2*(2.0*kt4+kt2*ktso2[mu][nu])+c12*(kt42+kt2*kt4*ktso2[mu][nu]+kt22*ktpo2[mu][nu]));
      A[mu][nu]-=(1.0 - kron_delta[mu][nu]);
    }

  return A;
}

vector<double> compute_Z(valarray<double>& Int,valarray<double>& IntS);

// vector<double> compute_V(valarray<double>& Int);
// double compute_S(valarray<double>& Int);


#endif
