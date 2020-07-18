AC_INIT([template],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([bin])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(bin/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_RANLIB

AC_PROG_CC
AC_PROG_CXX
AC_LANG([C++])

#silent automake
AM_SILENT_RULES([yes])

AX_CXX_COMPILE_STDCXX_11(noext,mandatory)
AX_CXXFLAGS_WARN_ALL

dnl #add Eigen path, search it and enable
dnl AC_ARG_WITH(eigen,
dnl 	AS_HELP_STRING([--with-eigen],[Specify where to find Eigen/Dense]),
dnl 	with_eigen=${withval}
dnl 	CPPFLAGS="$CPPFLAGS -I$with_eigen",
dnl 	with_eigen=no)
dnl AC_MSG_RESULT([with Eigen... ${with_eigen}])
dnl AC_LANG_PUSH([C++])
dnl AC_CHECK_HEADERS([Eigen/Dense])
dnl AC_LANG_POP([C++])
dnl if test "$ac_cv_header_Eigen_Dense" == "no"
dnl then
dnl 	AC_MSG_ERROR(["Unable to find Eigen"])
dnl else
dnl 	AC_DEFINE_UNQUOTED([USE_EIGEN],[1],[Enable Eigen])
dnl fi
dnl AC_MSG_RESULT(enabling Eigen ... ${enable_eigen})

#check immintrin.h
AC_CHECK_HEADER(immintrin.h,
	simd_inst_set_default=avx,
	simd_inst_set_default=none
	AC_DEFINE_UNQUOTED([DISABLE_X86_INTRINSICS],1,[Disable intriniscs]))

#simd set to be used
AC_ARG_WITH(simd-inst-set,
	AS_HELP_STRING([--with-simd-inst-set=set],[Select the set of SIMD instruction (avx [default if possible], none [otherwise], mmx or avx512)]),
	with_simd_inst_set="${withval}",
	with_simd_inst_set=avx)
case "$with_simd_inst_set" in
     none) CPPFLAGS_SIMD=""
     	  SIMD_INST_SET=NONE;;
     mmx) CPPFLAGS_SIMD="-mmmx"
     	  SIMD_INST_SET=MMX;;
     avx) CPPFLAGS_SIMD="-mavx"
     	  SIMD_INST_SET=AVX;;
     avx512) CPPFLAGS_SIMD="-mavx512f"
     	  SIMD_INST_SET=AVX512;;
      *) AC_MSG_ERROR(["Unkwnown SIMD instruction set ${withval}"])
esac

#check that we asked simd type only if available
if test SIMD_INST_SET != NONE && test "$ac_cv_header_immintrin_h" == ""
then
	AC_MSG_ERROR(["Cannot use simd type if immintrin.h is not found"])
fi

if test "$CPPFLAGS_SIMD" != ""
then
	CPPFLAGS="$CPPFLAGS_SIMD -mfma $CPPFLAGS"
fi

AC_MSG_RESULT(with SIMD insttruction set ... ${with_simd_inst_set})
AC_DEFINE_UNQUOTED([SIMD_INST_SET],[$SIMD_INST_SET],[Instruction set to be used])

#real type to be used
AC_ARG_WITH(real-type,
	AS_HELP_STRING([--with-real-type=type],[Select the type used for real number: float, double (default) or quadruple]),
	with_real_type="${withval}",
	with_real_type=double)
case "$with_real_type" in
     float) REAL=float;;
     double) REAL=double
     	     SIMD_SUFF=d;;
     quad) REAL=__float128
     	   LIBS="$LIBS -lquadmath";;
      *) AC_MSG_ERROR(["Unkwnown real type ${withval}"])
esac

if test "REAL" == "__float128" -a "$SIMD_INST_SET" != "NONE"
then
	AC_MSG_ERROR(["Cannot vectorize quadruple precision"])
fi

AC_MSG_RESULT(with real type ... ${with_real_type})
AC_DEFINE_UNQUOTED([REAL],[$REAL],[Data type for real number])
AC_DEFINE_UNQUOTED([SIMD_SUFF],[$SIMD_SUFF],[Suffix for simd data type])

#openMP
AC_OPENMP
CFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
CPPFLAGS="$CPPFLAGS $OPENMP_CPPFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"

AC_CONFIG_FILES(Makefile bin/Makefile)

AC_OUTPUT
