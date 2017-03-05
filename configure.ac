#############################################################################
# Configuration for ctools - tools for science analysis of CTA data         #
# ------------------------------------------------------------------------- #
# Copyright (C) 2011-2016 Juergen Knoedlseder                               #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
# ------------------------------------------------------------------------- #
# Process this file with autoconf to produce a configure script.            #
# ------------------------------------------------------------------------- #
# ToDo:                                                                     #
# - add possibility to build as Framework on Darwin                         #
# - strip "-Wstrict-prototypes" for C++ compiling                           #
#############################################################################

#############################################################################
# Initialisation                                                            #
#############################################################################
AC_INIT([ctools], [1.3.0.dev1], [jurgen.knodlseder@irap.omp.eu], [ctools])
AC_CONFIG_SRCDIR([src/ctselect/ctselect.cpp])
AC_CONFIG_MACRO_DIR([m4])
#AM_INIT_AUTOMAKE([gnu no-dependencies foreign])
AM_INIT_AUTOMAKE([gnu -Wno-portability foreign])
#AM_INIT_AUTOMAKE([gnu -Wno-portability foreign serial-tests])
AC_CONFIG_HEADER([config.h])


#############################################################################
# Library versioning (current:revision:age)                                 #
# ------------------------------------------------------------------------- #
# Here are a set of rules to update your library version information:       #
# 1. Update the version information only immediately before a public        #
#    release of your software. More frequent updates are unnecessary,       #
#    and only guarantee that the current interface number gets larger       #
#    faster.                                                                #
# 2. If the library source code has changed at all since the last update,   #
#    then increment revision.                                               #
# 3. If any interfaces have been added, removed, or changed since the last  #
#    update, increment current, and set revision to 0.                      #
# 4. If any interfaces have been added since the last public release, then  #
#    increment age.                                                         #
# 5. If any interfaces have been removed or changed since the last public   #
#    release, then set age to 0.                                            #
#############################################################################
CTOOLS_LT_VERSION="3:0:0"
AC_SUBST(CTOOLS_LT_VERSION)


#############################################################################
# Set prefix                                                                #
# ------------------------------------------------------------------------- #
# By default, ctools get installed into "/usr/local/gamma".  The default    #
# can be overwritten using the --prefix option.                             #
#############################################################################
AC_PREFIX_DEFAULT([/usr/local/gamma])
if test "x$prefix" = xNONE; then
  ctools_prefix="/usr/local/gamma"
else
  ctools_prefix=${prefix}
fi
AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], "${ctools_prefix}", [Location where library is installed])


#############################################################################
# Check for programs that are required by automake                          #
#############################################################################
AC_PROG_MAKE_SET


#############################################################################
# Determine the build and host system type                                  #
#############################################################################
AC_CANONICAL_HOST


#############################################################################
# Check for --enable-coverage option
#############################################################################
AC_MSG_CHECKING(whether to build ctools for code coverage measure)
AC_ARG_ENABLE([coverage],
              AS_HELP_STRING([--enable-coverage],
                             [turn on code coverage measure [default=no]]),
              [enable_coverage="$enableval"],
              [enable_coverage="no"])
if test "x$enable_coverage" = "xyes"; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi


#############################################################################
# Check for valid C++ compiler.  This has to come before the check for      #
# header files, otherwise we won't detect if the C++ compiler is missing.   #
#                                                                           #
# A special kluge was added for Mac OS X >=10.7. As Apple switched to clang #
# we search from this version on for the clang compiler.                    #
#                                                                           #
# We also signal whether Mac OS X >=10.9 (Mavericks) and Mac OS X >=10.11   #
# (El Capitan).                                                             #
#############################################################################
# Set default Mac OS X version flags
osx_lion=no
osx_mavericks=no
osx_el_capitan=no

# Set default compiler lists
cxx_list="g++ c++ gpp aCC CC cxx cc++ cl.exe FCC KCC RCC xlC_r xlC clang++"
cc_list="gcc cc cl.exe clang"

# Check if we have Mac OS X >= 10.7 and >= 10.11
if test `uname -s` = "Darwin"; then
  version=`sw_vers -productVersion | sed 's/\(10\.[[0-9]]*\).*/\1/' | sed -e 's/\.//g'`
  if [ test $version -ge 107 ]; then
    cxx_list="clang++ g++ c++ gpp aCC CC cxx cc++ cl.exe FCC KCC RCC xlC_r xlC"
    cc_list="clang gcc cc cl.exe"
    osx_lion=yes
  fi
  if [ test $version -ge 109 ]; then
    osx_mavericks=yes
  fi
  if [ test $version -ge 1011 ]; then
    osx_el_capitan=yes
  fi
fi
AC_MSG_CHECKING([for Mac OS X >= 10.7])
if [test "x$osx_lion" = "xyes"]; then
        AC_MSG_RESULT(yes)
else
        AC_MSG_RESULT(no)
fi
AC_MSG_CHECKING([for Mac OS X >= 10.9])
if [test "x$osx_mavericks" = "xyes"]; then
        AC_MSG_RESULT(yes)
else
        AC_MSG_RESULT(no)
fi
AC_MSG_CHECKING([for Mac OS X >= 10.11])
if [test "x$osx_el_capitan" = "xyes"]; then
        AC_MSG_RESULT(yes)
else
        AC_MSG_RESULT(no)
fi
AM_CONDITIONAL(IS_ELCAPITAN, test "x$osx_el_capitan" = "xyes")

# If coverage is selected then save current C++ options
if test "x$enable_coverage" = "xyes"; then
  SAVED_CXXFLAGS=$CXXFLAGS
  SAVED_CFLAGS=$CFLAGS
fi

# Our language is C++
AC_LANG(C++)
AC_PROG_CXX($cxx_list)
AC_PROG_CC($cc_list)

# If coverage is selected then restore former C++ options
if test "x$enable_coverage" = "xyes"; then
  CXXFLAGS=$SAVED_CXXFLAGS
  CFLAGS=$SAVED_CFLAGS
fi


#############################################################################
# Darwin option: Determine architecture of universal library to be built    #
# ------------------------------------------------------------------------- #
# --with-univeral-archs=VALUE, where VALUE is one of:                       #
# 32-bit - ppc, i386                                                        #
# 64-bit - ppc64, x86_64                                                    #
# 3-way  - ppc, i386, x86_64 (default)                                      #
# intel  - i386, x86_64                                                     #
# all    - ppc, ppc64, i386, x86_64                                         #
#############################################################################
UNIVERSAL_ARCHS="3-way"
AC_SUBST(LIPO_32BIT_FLAGS)
AC_MSG_CHECKING(for --with-universal-archs)
AC_ARG_WITH(universal-archs,
    AS_HELP_STRING([--with-universal-archs=ARCH], [select architectures for universal build ("32-bit", "64-bit", "3-way", "intel" or "all")]),
[
	AC_MSG_RESULT($withval)
	UNIVERSAL_ARCHS="$withval"
],
[
 	AC_MSG_RESULT(3-way)
])


#############################################################################
# Darwin option: Determine if a universal library should be built           #
# ------------------------------------------------------------------------- #
# --enable-universalsdk=PATH takes an optional argument that specifies      #
# which OSX SDK should be used to perform the build.                        #
# This defaults to:                                                         #
# "/Developer/SDKs/MacOSX.10.4u.sdk" for a '32-bit' build                   #
# "/Developer/SDKs/MacOSX.10.5.sdk"  for any other build                    #
#############################################################################
AC_MSG_CHECKING([for --enable-universalsdk])
AC_ARG_ENABLE(universalsdk,
	AS_HELP_STRING([--enable-universalsdk@<:@=SDKDIR@:>@], [Build against Mac OS X 10.4u SDK (ppc/i386)]),
[
	case $enableval in
	yes)
        if [test "${UNIVERSAL_ARCHS}" = "32-bit";] then
            enableval=/Developer/SDKs/MacOSX10.4u.sdk
        else
            enableval=/Developer/SDKs/MacOSX10.5.sdk
        fi
		if [test ! -d "${enableval}";] then
			enableval=/
        fi
		;;
	esac
	case $enableval in
	no)
		UNIVERSALSDK=
		enable_universalsdk=
		;;
	*)
		UNIVERSALSDK=$enableval
		if [test ! -d "${UNIVERSALSDK}"]; then
			AC_MSG_ERROR([--enable-universalsdk specifies non-existing SDK: ${UNIVERSALSDK}])
		fi
		;;
	esac
	
],[
   	UNIVERSALSDK=
	enable_universalsdk=
])
if [test -n "${UNIVERSALSDK}"]; then
	AC_MSG_RESULT(${UNIVERSALSDK})
else
	AC_MSG_RESULT(no)
fi
AC_SUBST(UNIVERSALSDK)


#############################################################################
# Darwin option: Check whether --with-universal-archs was specified without #
#                --enable-universalsdk                                      #
#############################################################################
if [test "${withval}"]; then
    if [test "${enable_universalsdk}"]; then
		:
	else
		AC_MSG_ERROR([--with-universal-archs without --enable-universalsdk.])
	fi
fi


#############################################################################
# Record MACOSX_DEPLOYMENT_TARGET                                           #
# ------------------------------------------------------------------------- #
# Record the configure-time value of MACOSX_DEPLOYMENT_TARGET, it may       #
# influence the way we can build extensions, so distutils needs to check it #
#############################################################################
AC_SUBST(CONFIGURE_MACOSX_DEPLOYMENT_TARGET)
AC_SUBST(EXPORT_MACOSX_DEPLOYMENT_TARGET)
CONFIGURE_MACOSX_DEPLOYMENT_TARGET=
EXPORT_MACOSX_DEPLOYMENT_TARGET='#'


#############################################################################
# Set name for machine-dependent library files                              #
#############################################################################
AC_SUBST(MACHDEP)
AC_MSG_CHECKING(MACHDEP)
if test -z "$MACHDEP"; then
    ac_sys_system=`uname -s`
	if test "$ac_sys_system" = "AIX" \
	-o "$ac_sys_system" = "UnixWare" -o "$ac_sys_system" = "OpenUNIX"; then
		ac_sys_release=`uname -v`
	else
		ac_sys_release=`uname -r`
	fi
	ac_md_system=`echo $ac_sys_system |
			   tr -d '[/ ]' | tr '[[A-Z]]' '[[a-z]]'`
	ac_md_release=`echo $ac_sys_release |
			   tr -d '[/ ]' | sed 's/^[[A-Z]]\.//' | sed 's/\..*//'`

    # Set MACHDEP
    MACHDEP="$ac_md_system$ac_md_release"
    case $MACHDEP in
    cygwin*) MACHDEP="cygwin";;
    darwin*) MACHDEP="darwin";;
    atheos*) MACHDEP="atheos";;
    irix646) MACHDEP="irix6";;
    '')	MACHDEP="unknown";;
    esac
fi
AC_MSG_CHECKING(machine type as reported by uname -m)
ac_sys_machine=`uname -m`
AC_MSG_RESULT($ac_sys_machine)


#############################################################################
# Set system-dependent library definitions                                  #
#############################################################################
AC_MSG_CHECKING([for the shared library paths environment variable])
case $ac_sys_system in
    BeOS*)
        SHAREDSUFFIX='.so'
        ;;
    CYGWIN*)
        SHAREDSUFFIX='.dll.a'
        ;;
    SunOS*)
        SHAREDSUFFIX='.so'
        RUNSHARED='LD_LIBRARY_PATH'
        ;;
    Linux*|GNU*|NetBSD*|FreeBSD*|DragonFly*)
        SHAREDSUFFIX='.so'
        RUNSHARED='LD_LIBRARY_PATH'
        ;;
    hp*|HP*)
        case `uname -m` in
            ia64)
                SHAREDSUFFIX='.so'
                ;;
            *)
                SHAREDSUFFIX='.sl'
                ;;
        esac
        RUNSHARED='SHLIB_PATH'
        ;;
    OSF*)
        SHAREDSUFFIX='.so'
        RUNSHARED='LD_LIBRARY_PATH'
        ;;
    atheos*)
        SHAREDSUFFIX='.so'
        RUNSHARED='DLL_PATH'
	    ;;
    Darwin*)
        SHAREDSUFFIX='.dylib'
        RUNSHARED='DYLD_LIBRARY_PATH'
        ;;
    AIX*)
        SHAREDSUFFIX='.so'
        RUNSHARED='LIBPATH'
        ;;
esac
AC_MSG_RESULT($RUNSHARED)
AC_SUBST(SHAREDSUFFIX)
AC_SUBST(RUNSHARED)


#############################################################################
# Set -arch flags for universal builds on OSX and tweak BASECFLAGS based on #
# compiler and platform (note that BASECFLAGS are not yet used).            #
#############################################################################
# Initialise the -arch flags for universal builds on OSX
UNIVERSAL_ARCH_FLAGS=
AC_SUBST(UNIVERSAL_ARCH_FLAGS)

# Tweak BASECFLAGS based on compiler and platform
case $GCC in
yes)

    # if using gcc on alpha, use -mieee to get (near) full IEEE 754
    # support.  Without this, treatment of subnormals doesn't follow
    # the standard.
    case $ac_sys_machine in
        alpha*)
            BASECFLAGS="$BASECFLAGS -mieee"
            ;;
    esac

    case $ac_sys_system in
	SCO_SV*)
	    BASECFLAGS="$BASECFLAGS -m486 -DSCO5"
	    ;;
	# is there any other compiler on Darwin besides gcc?
	Darwin*)
	    # -Wno-long-double, -no-cpp-precomp, and -mno-fused-madd
	    # used to be here, but non-Apple gcc doesn't accept them.
	    if test "${CC}" = gcc; then
            AC_MSG_CHECKING(which compiler should be used)
		    case "${UNIVERSALSDK}" in
		    */MacOSX10.4u.sdk)
			    # Build using 10.4 SDK, force usage of gcc when the 
			    # compiler is gcc, otherwise the user will get very
			    # confusing error messages when building on OSX 10.6
			    CC=gcc-4.0
			    CXX=g++-4.0
			    CPP=cpp-4.0
                ac_ct_CC="$CC"
                ac_ct_CXX="$CXX"
                ac_cv_prog_ac_ct_CC="$ac_ct_CC"
                ac_cv_prog_ac_ct_CXX="$ac_ct_CXX"
			    ;;
		    esac
		    AC_MSG_RESULT($CXX)
        fi

	    # Calculate the right deployment target for this build.
	    cur_target=`sw_vers -productVersion | sed 's/\(10\.[[0-9]]*\).*/\1/'`
	    if test ${cur_target} '>' 10.2; then
		    cur_target=10.3
		    if test ${enable_universalsdk}; then
			    if test "${UNIVERSAL_ARCHS}" = "all"; then
				    # Ensure that the default platform for a 
				    # 4-way universal build is OSX 10.5, 
				    # that's the first OS release where 
				    # 4-way builds make sense.
				    cur_target='10.5'
			    elif test "${UNIVERSAL_ARCHS}" = "3-way"; then
				    cur_target='10.5'
			    elif test "${UNIVERSAL_ARCHS}" = "intel"; then
				    cur_target='10.5'
			    elif test "${UNIVERSAL_ARCHS}" = "64-bit"; then
				    cur_target='10.5'
			    elif test "${UNIVERSAL_ARCHS}" = "32-bit"; then
				    cur_target='10.4'
			    fi
		    else
			    if test `/usr/bin/arch` = "i386"; then
				    # On Intel macs default to a deployment
				    # target of 10.4, that's the first OSX
				    # release with Intel support.
				    cur_target="10.4"
			    fi
		    fi
	    fi
	    CONFIGURE_MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET-${cur_target}}
	    
        # Make sure that deployment target is >= SDK
        case "${UNIVERSALSDK}" in
        */MacOSX10.4u.sdk)
            if [test ${cur_target} '>' 10.4;] then
	            AC_MSG_ERROR([Invalid MACOSX_DEPLOYMENT_TARGET=$CONFIGURE_MACOSX_DEPLOYMENT_TARGET for SDK 10.4. Specify --enable-universalsdk=/Developer/SDKs/MacOSX${cur_target}.sdk])
            fi
            ;;
        */MacOSX10.5.sdk)
            if [test ${cur_target} '>' 10.5;] then
	            AC_MSG_ERROR([Invalid MACOSX_DEPLOYMENT_TARGET=$CONFIGURE_MACOSX_DEPLOYMENT_TARGET for SDK 10.5. Specify --enable-universalsdk=/Developer/SDKs/MacOSX${cur_target}.sdk])
            fi
            ;;
        */MacOSX10.6.sdk)
            if [test ${cur_target} '>' 10.6;] then
	            AC_MSG_ERROR([Invalid MACOSX_DEPLOYMENT_TARGET=$CONFIGURE_MACOSX_DEPLOYMENT_TARGET for SDK 10.6. Specify --enable-universalsdk=/Developer/SDKs/MacOSX${cur_target}.sdk])
            fi
            ;;
        esac
	    
        # Set compile, link and pre-processor flags dependent on architecture
	    if test "${enable_universalsdk}"; then

	        # Make sure that MACOSX_DEPLOYMENT_TARGET is set in the 
	        # environment with a value that is the same as what we'll use
	        # in the Makefile to ensure that we'll get the same compiler
	        # environment during configure and build time.
	        MACOSX_DEPLOYMENT_TARGET="$CONFIGURE_MACOSX_DEPLOYMENT_TARGET"
	        export MACOSX_DEPLOYMENT_TARGET
	        EXPORT_MACOSX_DEPLOYMENT_TARGET=''

		    UNIVERSAL_ARCH_FLAGS=""
	        if test "$UNIVERSAL_ARCHS" = "32-bit" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch ppc -arch i386"
		        ARCH_RUN_32BIT=""
		        LIPO_32BIT_FLAGS=""
            elif test "$UNIVERSAL_ARCHS" = "64-bit" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch ppc64 -arch x86_64"
		        LIPO_32BIT_FLAGS=""
		        ARCH_RUN_32BIT="true"
            elif test "$UNIVERSAL_ARCHS" = "all" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch i386 -arch ppc -arch ppc64 -arch x86_64"
		        LIPO_32BIT_FLAGS="-extract ppc7400 -extract i386"
		        ARCH_RUN_32BIT="/usr/bin/arch -i386 -ppc"
            elif test "$UNIVERSAL_ARCHS" = "intel" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch i386 -arch x86_64"
		        LIPO_32BIT_FLAGS="-extract i386"
		        ARCH_RUN_32BIT="/usr/bin/arch -i386"
            elif test "$UNIVERSAL_ARCHS" = "3-way" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch i386 -arch ppc -arch x86_64"
		        LIPO_32BIT_FLAGS="-extract ppc7400 -extract i386"
		        ARCH_RUN_32BIT="/usr/bin/arch -i386 -ppc"
		    else
	            AC_MSG_ERROR([proper usage is --with-universal-arch=32-bit|64-bit|all|intel|3-way])
		    fi
		    CXXFLAGS="${UNIVERSAL_ARCH_FLAGS} ${CXXFLAGS}"
		    if test "${UNIVERSALSDK}" != "/"; then
			    CPPFLAGS="-isysroot ${UNIVERSALSDK} ${CPPFLAGS}"
			    LDFLAGS="-isysroot ${UNIVERSALSDK} ${LDFLAGS}"
			    CXXFLAGS="-isysroot ${UNIVERSALSDK} ${CXXFLAGS}"
		    fi
	    fi
	    ;;
	OSF*)
	    BASECFLAGS="$BASECFLAGS -mieee"
	    ;;
    esac
    ;;
*)
    case $ac_sys_system in
    OpenUNIX*|UnixWare*)
	BASECFLAGS="$BASECFLAGS -K pentium,host,inline,loop_unroll,alloca "
	;;
    OSF*)
	BASECFLAGS="$BASECFLAGS -ieee -std"
    	;;
    SCO_SV*)
	BASECFLAGS="$BASECFLAGS -belf -Ki486 -DSCO5"
	;;
    esac
    ;;
esac


#############################################################################
# Checks for debugging and profiling options                                #
#############################################################################
# Check for debugging option
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

# Check for profiling option
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


#############################################################################
# Set linker options                                                        #
#############################################################################
case $ac_sys_system/$ac_sys_release in
Darwin/*)
    # Use -undefined dynamic_lookup whenever possible (10.3 and later).
    # This allows an extension to be used in any Python
    #if test ${MACOSX_DEPLOYMENT_TARGET} '>' 10.2; then
        if test "${enable_universalsdk}"; then
            LDFLAGS="${UNIVERSAL_ARCH_FLAGS} -isysroot ${UNIVERSALSDK} ${LDFLAGS}"
        fi
    #fi
    ;;
esac


#############################################################################
# Checks for header files                                                   #
#############################################################################
AC_HEADER_STDC


#############################################################################
# Checks for installer                                                      #
#############################################################################
AC_PROG_INSTALL


#############################################################################
# Checks for libtool.  We need to do this after the compiler has been       #
# determined.                                                               #
#############################################################################
AC_PROG_LIBTOOL


#############################################################################
# Checks for typedefs, structures, and compiler characteristics             #
#############################################################################
AC_TYPE_SIZE_T
# AC_TYPE_LONG_LONG_INT
AC_C_CONST
AC_C_INLINE



#############################################################################
# Checks if Python bindings are to be build                                 #
# ------------------------------------------------------------------------- #
# Python binding requires Python and either the files tools_wrap.cpp and    #
# tools.py or swig to create these files.  The following conditionals are   #
# set:                                                                      #
# WITH_PYTHON   - Python bindings will be built                             #
# HAVE_SWIG     - swig is available                                         #
# HAVE_WRAPPERS - wrappers are available                                    #
#############################################################################
# Initialise flags
enable_python="no"
has_python="no"
has_python_h="no"
has_swig="no"
has_wrappers="no"
build_wrappers="no"

# Python binding enabled?
AC_ARG_ENABLE([python-binding],
              AS_HELP_STRING([--enable-python-binding],
                             [build python binding [default=yes]]),
              [ac_enable_python_binding="$enableval"],
              [ac_enable_python_binding="yes"])
if test "x$ac_enable_python_binding" = "xno"; then
  AC_MSG_NOTICE([Python binding for ctools will not be built])
  AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
else

  # Do we have Python and Python.h?
  AM_PATH_PYTHON([2.3],
                 [has_python="yes"],
                 [AC_MSG_WARN([Python >= 2.3 not found. Python >= 2.3 is required to build python binding. Python can be obtained from http://www.python.org])])

  # It would really be weired if Python has not been found, as barely every system has Python. Maybe
  # it was not found because the autoconf is so old, and in an early version of the AM_PATH_PYTHON
  # there were no [Action-if-found] and [Action-if-not-found] arguments, hence $has_python is possibly
  # badly set. However, the Python executable and version should be set correctly, and we can test
  # for this to see whether Python is really there or not.
  if test "x$has_python" = "xno"; then
    AM_PYTHON_CHECK_VERSION([$PYTHON], [${PYTHON_VERSION}],
                            [has_python="yes"
                             AC_MSG_NOTICE([Wow! You're running a pretty old system. Found Python nevertheless:-)])],
                            [AC_MSG_WARN([Python not found.])])
  fi

  # If we have Python then verify the presence of the header and set some flags
  if test "x$has_python" = "xyes"; then

    # Set opening and closing symbols to make setup.py compatible for versions 2 and 3
    if test ${PYTHON_VERSION} '<' 3; then
      PYTHON_RELATIVE_IMPORT=""
      PYTHON_FROM_DOT_=""
    else
      PYTHON_RELATIVE_IMPORT="."
      PYTHON_FROM_DOT_="from . "
    fi
    AC_SUBST(PYTHON_RELATIVE_IMPORT)  
    AC_SUBST(PYTHON_FROM_DOT_)  

    # Do we have distutils?
    AC_MSG_CHECKING([for the distutils Python package])
    ac_distutils_result=`$PYTHON -c "import distutils" 2>&1`
    if test -z "$ac_distutils_result"; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
    fi

    # If we have distutils then get python path from distutils ...
    if test -z "$ac_distutils_result"; then
      python_path=`$PYTHON -c "import distutils.sysconfig; \
                               print(distutils.sysconfig.get_python_inc());"`
      if test -n "${python_path}"; then
        AC_CHECK_HEADERS($python_path/Python.h, [has_python_h="yes"], [has_python_h="no"]) 
        python_path="-I$python_path"
      fi
      PYTHON_CPPFLAGS=$python_path
      AC_SUBST([PYTHON_CPPFLAGS])

      # If we have Mavericks then get the MACOSX_DEPLOYMENT_TARGET from Python
      # If the target is older then Mavericks then set the MACOSX_DEPLOYMENT_TARGET
      # to this value as otherwise we get trouble with name mangling (that changed)
      if [test "x$osx_mavericks" = "xyes"]; then
        AC_MSG_CHECKING([Python MACOSX_DEPLOYMENT_TARGET])
        python_target=`$PYTHON -c "from distutils import sysconfig; \
                                   print(sysconfig.get_config_vars()[['MACOSX_DEPLOYMENT_TARGET']]);" 2> /dev/null`
        python_compact_target=`echo $python_target | sed 's/\(10\.[[0-9]]*\).*/\1/' | sed -e 's/\.//g'`
        case $python_compact_target in
          (*[!0-9]*|'')
            AC_MSG_RESULT(unknown)
            ;;
          (*)
            AC_MSG_RESULT(${python_target})
            if [ test $python_compact_target -le 108 ]; then
              CONFIGURE_MACOSX_DEPLOYMENT_TARGET="$python_target"
              MACOSX_DEPLOYMENT_TARGET="$python_target"
              export MACOSX_DEPLOYMENT_TARGET
              EXPORT_MACOSX_DEPLOYMENT_TARGET=''
              LIBS="${LIBS} -lstdc++"
              AC_MSG_NOTICE([set MACOSX_DEPLOYMENT_TARGET to ${python_target}])
            fi
            ;;
        esac
      fi

    # ... otherwise check for includes relative to exectuable (obsolete?)
    else                             
                             
      # Determine include path for python
      AC_PATH_PROG(PY_INCLUDE, python)  
      PY_INCLUDE=`echo $PY_INCLUDE | sed -e "s/bin/include/"`
      PY_INCLUDE="$PY_INCLUDE$PYTHON_VERSION"
  
      # Check if we have Python.h
      AC_CHECK_HEADERS($PY_INCLUDE/Python.h, [has_python_h="yes"], [has_python_h="no"])
  
      # Signal any problem
      if test "x$has_python_h" = "xno"; then
        AC_MSG_WARN([No Python header file "$PY_INCLUDE/Python.h" found. Please install appropriate Python development package.])
      fi
    fi

  fi

  # Do we have swig?
  AC_PATH_PROG([SWIG], [swig])
  if test "$SWIG"; then
    has_swig="yes"
  fi
 
  # Do we have wrappers (tools_wrap.cpp and tools.py)?
  has_wrappers="yes"
  AC_CHECK_FILES([$srcdir/pyext/ctools/tools_wrap.cpp
                  $srcdir/pyext/ctools/tools.py],,
                 [has_wrappers="no"])
  if test "x$has_wrappers" = "xno"; then
    AC_MSG_WARN([Python wrapper missing. Requires swig for wrapper generation.])
  fi

  # Don't we have wrappers but we have swig
  if [test "x$has_wrappers" = "xno" -a "x$has_swig" = "xyes"]; then
    build_wrappers="yes"
  fi

  # Signal if Python is enabled
  if test "x$ac_enable_python_binding" = "xyes" -a "x$has_python" = "xyes" -a "x$has_python_h" = "xyes"; then
    if test "x$has_swig" = "xyes" -o "x$has_wrappers" = "xyes"; then
      enable_python="yes"
    fi
  fi

fi
AM_CONDITIONAL(WITH_PYTHON, test "x$enable_python" = "xyes")
AM_CONDITIONAL(HAVE_SWIG, test "x$has_swig" = "xyes")
AM_CONDITIONAL(HAVE_WRAPPERS, test "x$has_wrappers" = "xyes")
AM_CONDITIONAL(BUILD_WRAPPERS, test "x$build_wrappers" = "xyes")

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
    AC_MSG_WARN([Neither pyext/ctools/tools_wrap.cpp and pyext/ctools/tools.py nor swig has not been found])
    AC_MSG_WARN([Python binding for ctools cannot be built])
    AC_MSG_WARN([Either use a ctools distribution with tools_wrap.cpp and tools.py shipped or install swig on your system])
    AC_MSG_WARN([swig can be obtained from http://www.swig.org/])
  fi
fi


#############################################################################
# Checks for coverage option                                                #
#############################################################################
# Initialise Python coverage with empty string
PYTHON_COVERAGE=
COVERAGE_CXXFLAGS=

# If coverage option has been selected then set some variables
if test "x$enable_coverage" = "xyes"; then
  CXXFLAGS="--coverage $CXXFLAGS"
  COVERAGE_CXXFLAGS="--coverage"
  LDFLAGS="--coverage $LDFLAGS"
  if test "x$has_python" = "xyes"; then
    PYTHON_COVERAGE="coverage run --branch --append "
  fi
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi
AC_SUBST(PYTHON_COVERAGE)
AC_SUBST(COVERAGE_CXXFLAGS)


#############################################################################
# Checks for GammaLib library                                               #
#############################################################################
# Initialise GammaLib availability
has_gammalib="no"

# First check for GammaLib in the standard path using the AC_LINK_IFELSE macro.
# This allows also search in paths that are eventually specified by LDFAGS.
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
  for i in "$ctools_prefix" "$GAMMALIB" /opt/local /usr/local /usr ; do
    if test -f $i/lib64/libgamma.so -o -f $i/lib64/libgamma.dylib; then
      if test -d $i/include/gammalib; then
        LIBDIR_GAMMALIB=$i/lib64
        INCDIR_GAMMALIB=$i/include/gammalib
        break
      fi
    fi
    if test -f $i/lib/libgamma.so -o -f $i/lib/libgamma.dylib; then
      if test -d $i/include/gammalib; then
        LIBDIR_GAMMALIB=$i/lib
        INCDIR_GAMMALIB=$i/include/gammalib
        break
      fi
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
else
  AC_MSG_ERROR([No suitable GammaLib library found.
ctools require GammaLib for compilation, so please install GammaLib before or,
if you think that it is installed, set the GAMMALIB environment variable to the
root directory of GammaLib.])
fi


#############################################################################
# Checks for OpenMP checking option (G_OPENMP)                              #
#############################################################################
# Initialise OpenMP availability
openmp_ok=no
gammalib_with_openmp=no
darwin_kluge=no

# Detect if GammaLib has been compiled with OpenMP
AC_MSG_CHECKING(whether GammaLib has been compiled with OpenMP)
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <GObservations.hpp>]],,)],			     
               [gammalib_with_openmp=no;AC_MSG_RESULT(no)],
               [gammalib_with_openmp=yes;AC_MSG_RESULT(yes)])

# Check if openMP is available
m4_ifdef([AC_OPENMP],[AC_OPENMP], [AX_OPENMP([openmp_ok=yes],[openmp_ok=no])])
if test "x$ac_cv_prog_cxx_openmp" != "xunsupported" && test "x$ac_cv_prog_cxx_openmp" != "x"; then
  openmp_ok=yes
fi

# If OpenMP is available
if test "x$openmp_ok" = "xyes"; then

  # Special Linux test to detect shared libgomp libraries that can not be loaded
  # dynamically
  if test "$OPENMP_CXXFLAGS" = "-fopenmp" -a "$ac_sys_system" = "Linux"; then
    AC_MSG_CHECKING(whether we can dynamically load libgomp)
    SAVED_LIBS=$LIBS
    LIBS="-ldl"
    AC_TRY_RUN(AC_LANG_PROGRAM([#include <dlfcn.h>],
                               [if (dlopen("libgomp.so.1", RTLD_NOW) == 0) return 1;]),
               [AC_MSG_RESULT(yes)],
               [enable_openmp="no"
                AC_MSG_RESULT(no)
                AC_MSG_WARN([OpenMP detected but shared libgomp.so.1 library cannot be dlopen()ed.])
                AC_MSG_WARN([If you use a gcc 4.2.x compiler, this is a known feature. Try to use gcc >= 4.3.])])
    LIBS=$SAVED_LIBS
  fi

  # Special Darwin test to see whether OpenMP compiles without adding the
  # gomp_thread_attr global or whether we need to add this global. 
  if test "$ac_sys_system" = "Darwin"; then
    AC_MSG_CHECKING(whether we need to add gomp_thread_attr)
    SAVED_CXXFLAGS=$CXXFLAGS
    CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
    AC_TRY_RUN(AC_LANG_PROGRAM([#include <omp.h>],
                               [omp_set_num_threads(1);]),
               [AC_MSG_RESULT(no)],
               [AC_TRY_RUN(AC_LANG_PROGRAM([#include <omp.h>;
                                            #include <pthread.h>;
                                            pthread_attr_t gomp_thread_attr;],
                                           [omp_set_num_threads(1);]),
                           [AC_MSG_RESULT(yes)
                            darwin_kluge="yes"
                            AC_DEFINE([HAVE_OPENMP_DARWIN_KLUGE], [1], [Define if OpenMP Darwin kluge is needed])],
                           [AC_MSG_RESULT(no)
                            AC_MSG_WARN([Do not manage to compile OpenMP.])])])
    CXXFLAGS=$SAVED_CXXFLAGS
  fi

  # If gammalib is compiled with OpenMP, ctools must also be compiled with
  # openmp
  if test "x$gammalib_with_openmp" = "xyes"; then
    enable_openmp="yes"
  else
    # ... otherwise, we use it optionally
    AC_MSG_CHECKING(whether to enable OpenMP)
    AC_ARG_ENABLE([openmp],
                  AS_HELP_STRING([--enable-openmp],
                                 [turn on OpenMP support [default=yes]]),
                  [enable_openmp="$enableval"],
                  [enable_openmp="yes"])
  fi

fi

# Set compiler flags accordingly
if test "x$enable_openmp" = "xyes"; then
  CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
  AC_DEFINE([HAVE_OPENMP], [1], [Define if OpenMP is enabled])
  AC_MSG_RESULT(yes)
else
  enable_openmp=no
  AC_MSG_RESULT(no)
fi
AC_SUBST(OPENMP_CXXFLAGS)

# Set system specific settings
OPENMP_LIB_LDFLAGS=
#if test "x$enable_openmp" = "xyes"; then
#  case $ac_sys_system/$ac_sys_release in
#  Darwin/*)
#    OPENMP_LIB_LDFLAGS=-lgomp
#    ;;
#  esac
#fi
AC_SUBST(OPENMP_LIB_LDFLAGS)

# Signal if ctools cannot be compiled as gammalib requires OpenMP but ctools
# cannot be compiled with OpenMP
if test "x$enable_openmp" = "xno" && test "x$gammalib_with_openmp" = "xyes"; then
  AC_MSG_ERROR([The GammaLib library has been built with openMP support but OpenMP is not available on this system. Please install OpenMP or re-build GammaLib without OpenMP support.])
fi
AM_CONDITIONAL(OPENMP,test "x$enable_openmp" = "xyes")


#############################################################################
# Build static code?                                                        #
#############################################################################
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
AC_SUBST(LIB_LINK_OPTION)


#############################################################################
# Checks for Doxygen                                                        #
#############################################################################
DX_HTML_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_INIT_DOXYGEN([ctools], [doc/Doxyfile], [doc/doxygen])


#############################################################################
# Set list of output files to be created                                    #
#############################################################################
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/ctbin/Makefile
                 src/ctlike/Makefile
                 src/ctmodel/Makefile
                 src/ctobssim/Makefile
                 src/ctskymap/Makefile
                 src/ctselect/Makefile
                 src/ctexpcube/Makefile
                 src/ctpsfcube/Makefile
                 src/ctedispcube/Makefile
                 src/ctbkgcube/Makefile
                 src/ctmapcube/Makefile
                 src/ctcubemask/Makefile
                 src/cttsmap/Makefile
                 src/ctbutterfly/Makefile
                 src/ctulimit/Makefile
                 src/cterror/Makefile
                 src/support/Makefile
                 src/ctools-setup
                 test/Makefile
                 test/test_python_ctools.sh
                 test/test_python_cscripts.sh
                 test/test_examples.py
                 test/testing.py
                 pyext/Makefile
                 pyext/setup.py
                 pyext/ctools/__init__.py
                 cscripts/__init__.py
                 cscripts/Makefile
                 caldb/Makefile
                 models/Makefile
                 examples/Makefile
                 ctools.pc
                 testconf.sh],
                [if test -f test/test_python_ctools.sh; then chmod +x test/test_python_ctools.sh; fi;
                 if test -f test/test_python_cscripts.sh; then chmod +x test/test_python_cscripts.sh; fi;
                 if test -f test/test_examples.py; then chmod +x test/test_examples.py; fi;
                 if test -f testconf.sh; then chmod +x testconf.sh; fi])
AC_OUTPUT


#############################################################################
# Print configuration summary to console                                    #
#############################################################################
echo
echo "  ctools configuration summary"
echo "  =============================="
if test "x$has_gammalib" = "xyes"; then
  echo "  * GammaLib support             (yes)   ${LIBDIR_GAMMALIB} ${INCDIR_GAMMALIB}"
else
  echo "  - GammaLib support             (no)    no gammalib library found"
fi

# Dump Python binding information
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "x$has_python" = "xyes"; then
    if test "x$has_python_h" = "xyes"; then
      if test "x$has_wrappers" = "xyes"; then
        if test "x$has_swig" = "xyes"; then
          echo "  * Make Python binding          (yes)   use swig for updates"
        else
          echo "  * Make Python binding          (yes)   use wrappers"
        fi
      else
        if test "x$has_swig" = "xyes"; then
          echo "  * Make Python binding          (yes)   use swig for building"
        else
          echo "  - Make Python binding          (no)    swig required to build wrappers"
        fi
      fi
    else
      echo "  - Make Python binding          (no)    Python.h not installed"
    fi
  else
    echo "  - Make Python binding          (no)    Python not installed"
  fi
  if test "x$has_python" = "xyes"; then
    echo "  * Python                       (yes)"
  else
    echo "  - Python                       (no)"
  fi
  if test "x$has_python_h" = "xyes"; then
    echo "  * Python.h                     (yes)"
  else
    echo "  - Python.h                     (no)"
  fi
  if test "x$has_swig" = "xyes"; then
    echo "  * swig                         (yes)"
  else
    echo "  - swig                         (no)"
  fi
  if test "x$has_wrappers" = "xyes"; then
    echo "  * Python wrappers              (yes)"
  else
    echo "  - Python wrappers              (no)"
  fi
else
  echo "  - Make Python binding          (no)"
fi

# Dump Doxygen information
if test "$DX_DOXYGEN" != ""; then
  echo "  * Doxygen                      (yes)   $DX_DOXYGEN"
else
  echo "  - Doxygen                      (no)"
fi

# Dump compile options
if test "x$enable_openmp" = "xyes"; then
  echo "  * Enable OpenMP                (yes)   (default)"
else
  echo "  - Enable OpenMP                (no)"
fi
if test "x$enable_debug" = "xyes"; then
  echo "  * Compile in debug code        (yes)"
else
  echo "  - Compile in debug code        (no)    (default)"
fi
if test "x$enable_profiling" = "xyes"; then
  echo "  * Enable code for profiling    (yes)"
else
  echo "  - Enable code for profiling    (no)    (default)"
fi

# Dump next step information
echo
if test "x$has_gammalib" = "xyes"; then
  echo "Now type 'make'"
  echo
else
  AC_MSG_ERROR([No suitable GammaLib library found. ctools require GammaLib
for compilation, so please install GammaLib before or provide the
location of the library and header files using the LDFLAGS and CPPFLAGS
environnment variables.])
fi

