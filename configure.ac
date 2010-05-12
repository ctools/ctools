# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.64])
AC_INIT([ctatools], [0.1.0], [knodlseder@cesr.fr])
AC_CONFIG_SRCDIR([src/ctselect/ctselect.cpp])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL

# Checks for libraries.
if test "x$prefix" = xNONE; then
  export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:/usr/local/lib/pkgconfig"
else
  export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$prefix/lib/pkgconfig"
fi
PKG_CHECK_MODULES(GammaLib, gammalib >= 0.0,,
AC_MSG_ERROR(
[GammaLib library not found. Consider adjusting the PKG_CONFIG_PATH 
environment variable if you installed GammaLib in a non-standard prefix.]))
CFLAGS="$CFLAGS $GammaLib_CFLAGS"
CPPFLAGS="$CPPFLAGS $GammaLib_CFLAGS"
LDFLAGS="$LDFLAGS $GammaLib_LIBS"

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/ctselect/Makefile])
AC_OUTPUT
