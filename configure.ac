# Configuration for ctatools - a set of tools for science analysis of CTA data
#
# Copyright (C) 2011 Jurgen Knodlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Process this file with autoconf to produce a configure script

AC_INIT([ctatools], [0.2.3], [knodlseder@cesr.fr], [ctatools])
AC_CONFIG_SRCDIR([src/ctselect/ctselect.cpp])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE


dnl ==========
dnl Set prefix
dnl ==========
AC_PREFIX_DEFAULT([/usr/local/gamma])
if test "x$prefix" = xNONE; then
  ctatools_prefix="/usr/local/gamma"
else
  ctatools_prefix=${prefix}
fi
AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], "${ctatools_prefix}", [Location where library is installed])


dnl ==============================================
dnl Add pkgconfig information of install directory
dnl ==============================================
# PKG_CONFIG_PATH="${ctatools_prefix}/lib/pkgconfig:/opt/local/lib/pkgconfig:/usr/local/lib/pkgconfig:/usr/lib/pkgconfig:$PKG_CONFIG_PATH"
# export PKG_CONFIG_PATH


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
                             [turn on debugging [default=no]]),
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
                             [turn on profiling [default=no]]),
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


dnl =========================================================================
dnl Checks if Python bindings are to be build. Python binding requires Python
dnl and either the files gammalib_wrap.cpp and gammalib.py or swig to create
dnl these files. The following conditionals are set:
dnl  WITH_PYTHON - Python bindings will be built
dnl  HAVE_SWIG    - swig is available
dnl =========================================================================

# Initialise flags
enable_python="no"
has_python="no"
has_swig="no"
has_wrappers="no"

# Python binding enabled?
AC_ARG_ENABLE([python-binding],
              AS_HELP_STRING([--enable-python-binding],
                             [build python binding [default=yes]]),
              [ac_enable_python_binding="$enableval"],
              [ac_enable_python_binding="yes"])
if test "x$ac_enable_python_binding" = "xno"; then
  AC_MSG_NOTICE([Python binding for gammalib will not be built])
  AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
fi

# Do we have Python and Python.h?
AM_PATH_PYTHON([2.2],
               [has_python="yes"],
               [AC_MSG_WARN([Python >= 2.2 not found. Python >= 2.2 is required to build python binding. Python can be obtained from http://www.python.org])])
if test "x$has_python" = "xyes"; then

  # Determine include path for python
  AC_PATH_PROG(PY_INCLUDE, python)  
  PY_INCLUDE=`echo $PY_INCLUDE | sed -e "s/bin/include/"`
  PY_INCLUDE="$PY_INCLUDE$PYTHON_VERSION"
  
  # Check if we have Python.h
  AC_CHECK_HEADERS($PY_INCLUDE/Python.h, [has_python="yes"], [has_python="no"])
  
  # Mac OS X kluge: if test before was unsuccessful then check for headers in /usr/include
  if test "x$has_python" = "xno"; then
    AC_CHECK_HEADERS("/usr/include/python$PYTHON_VERSION/Python.h", [has_python="yes"], [has_python="no"])
  fi

  # Signal any problem
  if test "x$has_python" = "xno"; then
    AC_MSG_WARN([No Python header file "$PY_INCLUDE/Python.h" found. Please install appropriate Python development package.])
  fi
fi

# Do we have swig?
AC_PATH_PROG([SWIG], [swig])
if test "$SWIG"; then
  has_swig="yes"
fi
AM_CONDITIONAL(HAVE_SWIG, test "x$has_swig" = "xyes")
 
# Do we have wrappers (gammalib_wrap.cpp and gammalib.py)?
AC_CHECK_FILES([pyext/gammalib_wrap.cpp pyext/gammalib.py],,)
if test "x$ac_cv_file_pyext_gammalib_wrap_cpp" = "xyes"; then
  if test "x$ac_cv_file_pyext_gammalib_py" = "xyes"; then
    has_wrappers="yes"
  fi
fi

# Signal if Python is enabled
if test "x$ac_enable_python_binding" = "xyes" -a "x$has_python" = "xyes"; then
  if test "x$has_swig" = "xyes" -o "x$has_wrappers" = "xyes"; then
    enable_python="yes"
  fi
fi
AM_CONDITIONAL(WITH_PYTHON, test "x$enable_python" = "xyes")

# Set Python environment
if test "x$enable_python" = "xyes"; then

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
  if test "x$has_python_wrap" = "xno" -a "x$has_swig" = "xno"; then
    AC_MSG_WARN([Neither pyext/ctatools_wrap.cpp and pyext/ctatools.py nor swig has not been found])
    AC_MSG_WARN([Python binding for ctatools cannot be built])
    AC_MSG_WARN([Either use a ctatools distribution with ctatools_wrap.cpp and ctatools.py shipped or install swig on your system])
    AC_MSG_WARN([swig can be obtained from http://www.swig.org/])
  fi
fi


dnl ===========================
dnl Checks for GammaLib library
dnl ===========================

# Initialise GammaLib availability
has_gammalib="no"

# First check for GammaLib in the standard path using the AC_LINK_IFELSE macro.
# This allows also search in paths that are eventually specified by LDFAGS.
AC_LANG(C++)
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lgamma"
AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <GammaLib.hpp>],
                                [GSkyDir dir])],
               [has_gammalib="yes"
                AC_MSG_NOTICE([GammaLib found in standard path])],
               [AC_MSG_WARN([GammaLib not found in any standard path])])
LDFLAGS=$SAVED_LDFLAGS  

# If GammaLib has not been found in any standard path then search now in a
# number of predefined locations.
if test "x$has_gammalib" = "xno"; then
  for i in "$ctatools_prefix" /opt/local /usr/local /usr ; do
    if test -f $i/lib64/libgamma.so -o -f $i/lib64/libgamma.so -a -d $i/include/gammalib; then
      LIBDIR_GAMMALIB=$i/lib64
      INCDIR_GAMMALIB=$i/include/gammalib
      break
    fi
    if test -f $i/lib/libgamma.so -o -f $i/lib/libgamma.dylib -a -d $i/include/gammalib; then
      LIBDIR_GAMMALIB=$i/lib
      INCDIR_GAMMALIB=$i/include/gammalib
      break
    fi
  done
  if test -n "$LIBDIR_GAMMALIB" -a -n "$INCDIR_GAMMALIB"; then
    has_gammalib="yes"
    LDFLAGS="$LDFLAGS -L$LIBDIR_GAMMALIB"
    CPPFLAGS="$CPPFLAGS -I$INCDIR_GAMMALIB"
    AC_MSG_NOTICE([GammaLib found in ${LIBDIR_GAMMALIB} and ${INCDIR_GAMMALIB}.])
  fi
fi

# Set GammaLib information
if test "x$has_gammalib" = "xyes"; then
  LIBS_GAMMALIB="-lgamma"
  LIBS="${LIBS} ${LIBS_GAMMALIB}"
  AC_DEFINE([HAVE_LIBGAMMA], [1], [Define if gammalib is available])
  AC_SUBST(LIBS_GAMMALIB)
  AC_SUBST(LIBDIR_GAMMALIB)
  AC_SUBST(INCDIR_GAMMALIB)
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
                 test/Makefile
                 pyext/Makefile
                 pyext/setup.py
                 ctatools.pc])
AC_OUTPUT


dnl ===========================
dnl Print configuration summary
dnl ===========================
echo
echo "  ctatools configuration summary"
echo "  =============================="
if test "x$has_gammalib" = "xyes"; then
  echo "  * GammaLib support             (yes)   ${LIBDIR_GAMMALIB} ${INCDIR_GAMMALIB}"
else
  echo "  - GammaLib support             (no)    no gammalib library found"
fi
if test ! "$SWIG"; then
  echo "  - swig                         (no)"
else
  echo "  * swig                         (yes)"
fi
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "x$has_python" = "xyes"; then
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


