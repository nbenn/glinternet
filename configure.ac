#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([glinternet], [0.9.1], [nbennett@ethz.ch])
AC_CONFIG_SRCDIR([src/c_routines.c])

# Find R home and set CC/CFLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
else
	echo "using R: ${R_HOME}"
fi
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`;
CFLAGS=`"${RBIN}" CMD config CFLAGS`
LIBS="${PKG_LIBS}"

OMPP=""
AC_ARG_WITH([OMPP], 
	AC_HELP_STRING(
  	[--with-OMPP=PATH],
  	[specify the path of kinst-ompp-papi]
  ),
  [OMPP=$withval]
)

WARNINGS=""
if test -n "${OMPP}"; then
	CC="$OMPP $CC"
else
	WARNINGS="-Wno-unknown-pragmas"
fi

echo "using CC: ${CC}"

# Checks for programs.
AC_PROG_CC
AC_OPENMP

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_RESTRICT
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt])

# since some systems have broken OMP libraries
# we also check that the actual package will work
ac_pkg_openmp=no
if test -n "${OPENMP_CFLAGS}"; then
  AC_MSG_CHECKING([whether OpenMP will work in a package])
  AC_LANG_CONFTEST(
  [AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
  PKG_CFLAGS="${OPENMP_CFLAGS}" PKG_LIBS="${OPENMP_CFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
  AC_MSG_RESULT([${ac_pkg_openmp}])
fi

# if ${ac_pkg_openmp} = "yes" then we have OMP, otherwise it will be "no"
if test "${ac_pkg_openmp}" = no; then
  OPENMP_CFLAGS=''
  # you could put AC_MSG_ERROR here is OpenMP is required
fi

GLNT_CFLAGS=""
GLNT_LIBS="-L/usr/lib64 -lnuma -lrt"

AC_SUBST(GLNT_CFLAGS)
AC_SUBST(GLNT_LIBS)
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(CC)
AC_SUBST(WARNINGS)

AC_CONFIG_FILES([src/Makefile])

AC_OUTPUT
