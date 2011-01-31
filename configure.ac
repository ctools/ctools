# Configuration for ctatools - a set of tools for science analysis of CTA data
#
# Copyright (C) 2011
#
# Process this file with autoconf to produce a configure script

AC_INIT([ctatools], [0.1.0], [knodlseder@cesr.fr], [ctatools], [])
AC_CONFIG_SRCDIR([src/ctselect/ctselect.cpp])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE


dnl ==========
dnl Set prefix
dnl ==========
AC_PREFIX_DEFAULT([/usr/local/gamma])
if test "x$prefix" = xNONE; then
  AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], ["/usr/local/gamma"],
                     [Location where library is installed])
else
  AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], "${prefix}",
                     [Location where library is installed])
fi


dnl ===========================
dnl things required by automake
dnl ===========================
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET


dnl =========================
dnl Check for system services
dnl =========================
AC_CANONICAL_HOST


dnl =======================
dnl Checks for header files
dnl =======================
AC_STDC_HEADERS


dnl =============================================================
dnl Checks for typedefs, structures, and compiler characteristics
dnl =============================================================
AC_TYPE_SIZE_T
AC_C_CONST
AC_C_INLINE


dnl ===================
dnl Checks for programs
dnl ===================
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL


dnl =================================
dnl Checks for debugging or profiling
dnl =================================
AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE([debug],
              AS_HELP_STRING([--enable-debug],
                             [turn on debugging (default=no)]),
              [enable_debug="$enableval"],
              [enable_debug="no"])
if test "x$enable_debug" = "xyes"; then
  CFLAGS="$CFLAGS -g"
  CXXFLAGS="$CXXFLAGS -g"
  FFLAGS="$FFLAGS -g"
  AC_DEFINE([G_DEBUG], [1], [Define if debugging is turned on])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi
AC_MSG_CHECKING(whether to enable profiling)
AC_ARG_ENABLE([profiling],
              AS_HELP_STRING([--enable-profiling],
                             [turn on profiling (default=no)]),
              [enable_profiling="$enableval"],
              [enable_profiling="no"])
if test "x$enable_profiling" = "xyes"; then
  CFLAGS="$CFLAGS -pg"
  CXXFLAGS="$CXXFLAGS -pg"
  FFLAGS="$FFLAGS -pg"
  AC_DEFINE([G_PROFILE], [1], [Define if profiling is turned on])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi


dnl ====================
dnl Checks for libraries
dnl ====================

# Get GammaLib configuration
if test "x$prefix" = xNONE; then
  export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:/usr/local/gamma/lib/pkgconfig"
else
  export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$prefix/lib/pkgconfig"
fi
PKG_CHECK_MODULES(GammaLib, gammalib >= 0.0,,
                  AC_MSG_ERROR([GammaLib library not found. Consider adjusting the PKG_CONFIG_PATH 
                                environment variable if you installed GammaLib in a non-standard prefix.]))
CFLAGS="$CFLAGS $GammaLib_CFLAGS"
CPPFLAGS="$CPPFLAGS $GammaLib_CFLAGS"
LDFLAGS="$LDFLAGS $GammaLib_LIBS"

# Check for GammaLib library and determine path to library
#AC_ARG_WITH([gammalib-lib-dir],
#            AC_HELP_STRING([--with-gammalib-lib-dir],
#                           [GammaLib library directory [[/usr/(local/(gamma/))lib(64)]]]),
#            [gammalib_lib_dir=$withval])
#if test "x$gammalib_lib_dir" != "x"; then
#  test -f $gammalib_lib_dir/libgamma.a -o -f $gammalib_lib_dir/libgamma.so -o -f $gammalib_lib_dir/libgamma.dylib && LIBDIR_GAMMALIB=$gammalib_lib_dir
#else
#  for i in /usr/local/gamma/lib /usr/local/gamma/lib64 /usr/local/lib /usr/local/lib64 /usr/lib /usr/lib64 ; do
#    test -f $i/libgamma.a -o -f $i/libgamma.so -o -f $i/libgamma.dylib && LIBDIR_GAMMALIB=$i
#  done
#fi
#if test -z "$LIBDIR_GAMMALIB"; then
#  if test "x$gammalib_lib_dir" != "x"; then
#    AC_MSG_ERROR(Cannot find GammaLib library under $gammalib_lib_dir)
#  else
#    AC_MSG_ERROR(Cannot find GammaLib library. Use --with-gammalib-lib-dir= to specify non-default path.)
#  fi
#fi
#AC_SUBST(LIBDIR_GAMMALIB)


dnl ========================================================================
dnl Checks for python binding. Python binding requires Python and either the
dnl files gammalib_wrap.cpp and gammalib.py or SWIG to create these files.
dnl The following conditionals are set:
dnl  HAVE_PYTHON           - Python >= 2.2 is available
dnl  HAVE_SWIG             - swig is available
dnl  ENABLE_PYTHON_BINDING - Python binding is requested
dnl  WITH_PYTHON_BINDING   - Either has swig or python wraper files
dnl ========================================================================
# Checks for Python (HAVE_PYTHON)
AM_PATH_PYTHON([2.2],
               [],
               [AC_MSG_WARN([Python not found. Python is required to build python binding. Python can be obtained from http://www.python.org])])
AM_CONDITIONAL(HAVE_PYTHON, test "$PYTHON")

# Checks for SWIG (HAVE_SWIG)
AC_PATH_PROG([SWIG], [swig])
AM_CONDITIONAL(HAVE_SWIG, test "$SWIG")

# Python binding enabled?
AC_ARG_ENABLE([python-binding],
              AS_HELP_STRING([--enable-python-binding],
                             [build python binding (default=yes)]),
              [ac_enable_python_binding="$enableval"],
              [ac_enable_python_binding="yes"])
AM_CONDITIONAL(ENABLE_PYTHON_BINDING, test "x$ac_enable_python_binding" = "xyes")

# Checks for ctatools_wrap.cpp and ctatools.py
AC_CHECK_FILES([pyext/ctatools_wrap.cpp pyext/ctatools.py],,)
has_python_wrap="no"
if test "x$ac_cv_file_pyext_ctatools_wrap_cpp" = "xyes"; then
  if test "x$ac_cv_file_pyext_ctatools_py" = "xyes"; then
    has_python_wrap="yes"
  fi
fi
AM_CONDITIONAL(WITH_PYTHON_BINDING, test "x$has_python_wrap" = "xyes" -o "$SWIG")

# Inform about missing components and set host dependent variables
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "$PYTHON"; then
    if test "x$has_python_wrap" != "xyes"; then
      if test ! "$SWIG"; then
        AC_MSG_WARN([Neither pyext/ctatools_wrap.cpp and pyext/ctatools.py nor swig has not been found])
        AC_MSG_WARN([Python binding for ctatools cannot be built])
        AC_MSG_WARN([Either use a ctatools distribution with ctatools_wrap.cpp and ctatools.py shipped or install swig on your system])
        AC_MSG_WARN([swig can be obtained from http://www.swig.org/])
      fi
    fi

    # Set host specific build prefix (needed for Solaris pycc bug)
    PYTHON_BUILD_PREFIX=
    case $host in
      sparc*-sun-solaris*)  PYTHON_BUILD_PREFIX="CC=CC" ;;
      *)                    PYTHON_BUILD_PREFIX="" ;;
    esac                    
    AC_SUBST(PYTHON_BUILD_PREFIX)
    
    # Set host specific extra libraries. Add all host specific libraries here that
    # need to be linked in the python module
    PYTHON_EXTRA_LIBS=
    case $host in
      sparc*-sun-solaris*)  PYTHON_EXTRA_LIBS="Cstd" ;;
      *)                    PYTHON_EXTRA_LIBS="" ;;
    esac                    
    AC_SUBST(PYTHON_EXTRA_LIBS)

  else
    AC_MSG_WARN([Python binding for ctatools cannot be built as Python seems not to be installed on your system])
  fi
else
  AC_MSG_NOTICE([Python binding for ctatools will not be built])
  AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
fi


dnl ==================
dnl Build static code?
dnl ==================
if test "x$enable_profiling" = "xyes"; then
  AC_MSG_NOTICE([Link executables statically.])
  LIB_LINK_OPTION="-all-static"
else
  AC_MSG_CHECKING(whether to link executables statically)
  AC_ARG_ENABLE([static_link],
                AS_HELP_STRING([--enable-static-link],
                               [link executables statically (default=no)]),
                [static_link="$enableval"],
                [static_link="no"])
  if test "x$static_link" = "xyes"; then
    AC_MSG_RESULT(yes)
    LIB_LINK_OPTION="-all-static"
  else
    AC_MSG_RESULT(no)
  fi
fi
AC_SUBST(LIB_LINK_OPTION)


dnl ==========================
dnl Create configuration files
dnl ==========================
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/ctselect/Makefile
                 src/ctbin/Makefile
                 src/ctlike/Makefile
                 src/ctobssim/Makefile
                 pyext/Makefile
                 pyext/setup.py])
AC_OUTPUT


dnl ===========================
dnl Print configuration summary
dnl ===========================
echo
echo "  ctatools configuration summary"
echo "  =============================="
if test ! "$SWIG"; then
  echo "  - swig                         (no)"
else
  echo "  * swig                         (yes)"
fi
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "$PYTHON"; then
    if test "x$has_python_wrap" = "xyes"; then
      if test ! "$SWIG"; then
        echo "  * Python binding               (yes)"
      else
        echo "  * Python binding               (yes)   use swig for updates"
      fi
    else
      if test ! "$SWIG"; then
        echo "  - Python binding               (no)    swig not installed"
      else
        echo "  * Python binding               (yes)   use swig for building"
      fi
    fi
  else
    echo "  - Python binding               (no)    python not installed"
  fi
else
  echo "  - Python binding               (no)"
fi
echo
echo "Now type 'make'"
echo
