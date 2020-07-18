#ifndef _REAL_HPP
#define _REAL_HPP

#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <quadmath.h>

using Real=REAL;

#if REAL == __float128

#define IMPORT_128_UNARY_FUNCTION(NAME)		\
  inline __float128  NAME (__float128 x)	\
  {						\
    return NAME ## q (x);			\
  }

IMPORT_128_UNARY_FUNCTION(fabs)
IMPORT_128_UNARY_FUNCTION(log)
IMPORT_128_UNARY_FUNCTION(sin)
IMPORT_128_UNARY_FUNCTION(cos)

#undef IMPORT_128_UNARY_FUNCTION

#endif

#endif
