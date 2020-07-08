#ifndef FUNC_HPP
#define FUNC_HPP


using namespace Eigen;
using namespace std;

inline Vector4d cwiseProduct(const Vector4d& a, const Vector4d& b)
{
  return a.cwiseProduct(b);
}
inline Vector4d cwiseProduct(const Vector4d& a, const Vector4d& b, const Vector4d& c)
{
  return a.cwiseProduct(b).cwiseProduct(c);
}
inline Vector4d cwiseProduct(const Vector4d& a, const Vector4d& b, const Vector4d& c, const Vector4d& d)
{
  return a.cwiseProduct(b).cwiseProduct(c).cwiseProduct(d);
}
inline Vector4d cwiseProduct(const Vector4d& a, const Vector4d& b, const Vector4d& c, const Vector4d& d, const Vector4d& e)
{
  return a.cwiseProduct(b).cwiseProduct(c).cwiseProduct(d).cwiseProduct(e);
}

inline double prod(const Vector4d& a)
{
  return a.dot(Id);
}
inline double prod(const Vector4d& a,const Vector4d& b)
{
  return a.dot(b);
};
inline double prod(const Vector4d& a,const Vector4d& b,const Vector4d& c)
{
  return cwiseProduct(a,b,c).dot(Id);
};
inline double prod(const Vector4d& a,const Vector4d& b,const Vector4d& c,const Vector4d& d)
{
  return cwiseProduct(a,b,c,d).dot(Id);
};
inline double prod(const Vector4d& a,const Vector4d& b,const Vector4d& c,const Vector4d& d,const Vector4d& e)
{
  return cwiseProduct(a,b,c,d,e).dot(Id);
};

inline Matrix4d make_propagator(const Vector4d& kt2_dir, const Vector4d& kt4_dir, const Vector4d& kt6_dir)
{
  Matrix4d ktpo2, ktso2, A;

  double c12 = c1*c1, c13 = c12*c1;
  int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

  double kt2 = prod(kt2_dir), kt4 = prod(kt4_dir), kt6 = prod(kt6_dir);

  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++){

      ktpo2(mu,nu)=1.0;
      ktso2(mu,nu)=0.0;

      for(int rho=0;rho<4;rho++){
        if(mu!=rho && nu!=rho)
        {
          ktpo2(mu,nu)*=kt2_dir(rho);
          ktso2(mu,nu)+=kt2_dir(rho);
        }
      }
    }

  double kt22=kt2*kt2;
  double kt23=kt2*kt2*kt2;
  double kt42=kt4*kt4;

  double Deltakt=(kt2-c1*kt4)*(kt2-c1*(kt22+kt4)+0.5*c12*(kt23+2*kt6-kt2*kt4));
  for(int rho=0;rho<4;rho++) Deltakt-=4*c13*kt4_dir(rho)*ktpo2(rho,rho);

  for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++){
      A(mu,nu)=(1.0 - kron_delta[mu][nu])/Deltakt*(kt22-c1*kt2*(2.0*kt4+kt2*ktso2(mu,nu))+c12*(kt42+kt2*kt4*ktso2(mu,nu)+kt22*ktpo2(mu,nu)));
      A(mu,nu)-=(1.0 - kron_delta[mu][nu]);
    }

  return A;
}

// vector<double> make_vector(vector<double>& np)
// {
//   for (int i=0; i<4; i++) np[i] *= 2.0*M_PI/V[i];
//   return np;
// }


#endif
