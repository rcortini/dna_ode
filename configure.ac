#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.66])
AC_INIT([dna_ode], [0.0.0])
AC_CONFIG_HEADERS([config.h])
PACKAGE=DNAODE
VERSION=0.0.0

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_STDC
AM_PROG_CC_C_O

# init automake and libtool
AM_INIT_AUTOMAKE(subdir-objects)
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4])

# Checks for libraries.

# 1. Checks for GSL
AC_ARG_WITH(gsl,
            [  --with-gsl=DIR        Directory where the GSL is installed (optional)],
            GSL_DIR="$withval")

GSL_CONFIG=
AC_PATH_PROG([GSL_CONFIG], [gsl-config],, [$GSL_DIR/bin:$PATH])

GSL_INCLUDES=
GSL_LIBS=

if test -x "$GSL_CONFIG"; then
  GSL_INCLUDES=`$GSL_CONFIG --cflags`
  GSL_LIBS=`$GSL_CONFIG --libs`
else
  AC_ERROR([GSL not found. You need gsl to compile this package.
      See http://www.gnu.org/software/gsl/])
fi

AC_SUBST(GSL_INCLUDES)
AC_SUBST(GSL_LIBS)

# 2. Check for ODE
PKG_CHECK_MODULES(ODE, ode)

# 3. Check for libconfig
PKG_CHECK_MODULES(CONFIG, libconfig)

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_LIB(m, [floor,sqrt,round,sin,cos,acos])

# Global CFLAGS
WARN_FLAGS="-Wall -Wextra -Wshadow -pedantic-errors"
CFLAGS="$WARN_FLAGS -O2"
CFLAGS_save="$CFLAGS"

# enable "optimize option", which will compile using only the -O3 option
AC_ARG_ENABLE([optimize],
    [AC_HELP_STRING([--enable-optimize],
	    [disable debugging (removes -g flags), disable compiler warnings\
	     and set optimization level O3 (disabled)])],
    [CFLAGS="-O3"], [CFLAGS="$CFLAGS_save"])

# enable "more-debug", which will compile using extra debugging options
CFLAGS_save="$CFLAGS"
AC_ARG_ENABLE([more-debug],
    [AC_HELP_STRING([--enable-more-debug],
	    [enable more debugging information (disabled)])],
     [CFLAGS="$WARN_FLAGS -g3 -ggdb -O0"], [CFLAGS="$CFLAGS_save"])

# Global LDFLAGS

# output Makefiles
AC_CONFIG_FILES([Makefile
                 dna_ode/Makefile
                 dna_ode/core/Makefile
                 dna_ode/fdna/Makefile
                 dna_ode/sdna/Makefile
                 examples/.examples_in/Makefile
                 examples/Makefile])

# finalisation
AC_OUTPUT
