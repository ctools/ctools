# Process this file with autoconf to produce a configure script.

AC_INIT([ctatools], [0.1.0], [knodlseder@cesr.fr])
AC_CONFIG_SRCDIR([src/ctselect/ctselect.cpp])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Set prefix
AC_PREFIX_DEFAULT([/usr/local/gamma])
AC_DEFINE([PREFIX], [$prefix], 'Installation prefix')

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL

# Checks for libraries.
if test "x$prefix" = xNONE; then
  export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:/usr/local/gamma/lib/pkgconfig"
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

# Checks for SWIG
AC_PATH_PROG([SWIG], [swig])
AM_CONDITIONAL(HAVE_SWIG, test "$SWIG")

# Checks for Python
AM_PATH_PYTHON([2.2],
	       [],
	       [AC_MSG_WARN([Python not found. Python is required to build python binding. Python can be obtained from http://www.python.org])])
AM_CONDITIONAL(HAVE_PYTHON, test "$PYTHON")

AC_ARG_ENABLE([python-binding],
	      AS_HELP_STRING([--enable-python-binding],
			     [build python binding (default=yes)]),
			     [ac_enable_python_binding=no],
			     [ac_enable_python_binding=yes])
if test "x$ac_enable_python_binding" = "xyes"
then
    if test ! "$SWIG" -a ! "$PYTHON"
    then
	AC_MSG_WARN([Python binding for ctatools cannot be built. Ensure that SWIG and Python packages are available.])
    fi
else
    AC_MSG_NOTICE([Python binding for ctatools will not be built.])
    AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
fi
AM_CONDITIONAL(ENABLE_PYTHON_BINDING, test "x$ac_enable_python_binding" = "xyes")

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/ctselect/Makefile
                 src/ctbin/Makefile
                 src/ctlike/Makefile
                 pyext/Makefile
                 pyext/setup.py])
AC_OUTPUT
