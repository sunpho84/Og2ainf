#ifndef _SIMD_HPP
#define _SIMD_HPP

#ifdef HAVE_CONFIG_H
#include <config.hpp>
#endif

#include <array>

#include <immintrin.h>

#include "preprocessor.hpp"
#include "real.hpp"

#define MMX 2
#define AVX 4
#define AVX512 8

#if SIMD_INST_SET == AVX512
 #define SIMD_NBITS 512
#elif SIMD_INST_SET == AVX
 #define SIMD_NBITS 256
#elif SIMD_INST_SET == MMX
 #define SIMD_NBITS 128
#endif

#if SIMD_INST_SET != NONE

using vtype=CONCAT3(__m,SIMD_NBITS,SIMD_SUFF);
constexpr int N=sizeof(vtype)/sizeof(Real);

#else

using vtype=std::array<Real,1>;
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
