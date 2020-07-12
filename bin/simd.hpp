#ifndef _SIMD_HPP
#define _SIMD_HPP

#include <array>

#include <immintrin.h>

// passabile da configure.ac, eventualmente
#define USE_AVX

#if defined USE_AVX512
using vtype=__m512d;
constexpr int N=8;
#elif defined USE_AVX
using vtype=__m256d;
constexpr int N=4;
#elif defined USE_MMX
using vtype=__m128d;
constexpr int N=2;
#else
using vtype=std::array<double,1>;
constexpr int N=1;

inline vtype& operator+=(vtype& l,const vtype& r)
{
  l[0]+=r[0];
  return l;
}

inline vtype& operator-=(vtype& l,const vtype& r)
{
  l[0]-=r[0];
  return l;
}

inline vtype& operator-=(vtype& l,const double& r)
{
  l[0]-=r;
  return l;
}

inline vtype& operator*=(vtype& l,const vtype& r)
{
  l[0]*=r[0];
  return l;
}

inline vtype operator*(const vtype& l,const vtype& r)
{
  vtype o;
  o[0]=l[0]*r[0];
  return o;
}

inline vtype operator*(const double& l,const vtype& r)
{
  vtype o;
  o[0]=l*r[0];
  return o;
}

inline vtype operator+(const vtype& l,const vtype& r)
{
  vtype o;
  o[0]=l[0]+r[0];
  return o;
}

inline vtype operator+(const double& l,const vtype& r)
{
  vtype o;
  o[0]=l+r[0];
  return o;
}

inline vtype operator-(const vtype& l,const vtype& r)
{
  vtype o;
  o[0]=l[0]-r[0];
  return o;
}

inline vtype operator-(const vtype& l,const double& r)
{
  vtype o;
  o[0]=l[0]-r;
  return o;
}

inline vtype operator-(const double& l,const vtype& r)
{
  vtype o;
  o[0]=l-r[0];
  return o;
}

inline vtype operator*(const vtype& l,const double& r)
{
  vtype o;
  o[0]=l[0]*r;
  return o;
}

inline vtype operator/(const vtype& l,const vtype& r)
{
  vtype o;
  o[0]=l[0]/r[0];
  return o;
}

inline vtype operator/(const double& l,const vtype& r)
{
  vtype o;
  o[0]=l/r[0];
  return o;
}

inline vtype operator/(const vtype& l,const double& r)
{
  vtype o;
  o[0]=l[0]/r;
  return o;
}

inline vtype operator-(const vtype& r)
{
  vtype o;
  o[0]=-r[0];
  return o;
}

#endif

#endif
