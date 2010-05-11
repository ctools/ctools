dnl =============================================================
dnl Process this file with autoconf to produce a configure script
dnl =============================================================
AC_INIT([CTATools], [0.1], [knodlseder@cesr.fr])
AC_CONFIG_SRCDIR(src/Makefile.am)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([gnu no-dependencies])
AM_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE


dnl ===========================
dnl things required by automake
dnl ===========================
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET


dnl ======================
dnl Check for which system
dnl ======================
dnl AC_CANONICAL_HOST


dnl ===================
dnl Checks for programs
dnl ===================
AC_PROG_CXX
AC_PROG_INSTALL
dnl AC_PROG_RANLIB
AC_PROG_LIBTOOL


dnl =======================
dnl Check compiler features
dnl =======================
dnl AC_TYPE_SIZE_T
dnl AC_C_CONST
dnl AC_C_VOLATILE
dnl AC_C_INLINE


dnl =======================
dnl Checks for header files
dnl =======================


dnl =============================================================
dnl Checks for typedefs, structures, and compiler characteristics
dnl =============================================================


dnl ============================
dnl Checks for library functions
dnl ============================


dnl ===================
dnl Checks for gammalib
dnl ===================
AC_SEARCH_LIBS(str, gammalib,
    AC_DEFINE([G_GAMMALIB], [1], [Define if gammalib is available]), ,[-lm])


dnl ======================================
dnl Set list of output files to be creates
dnl ======================================
AC_OUTPUT([Makefile
           src/Makefile
           src/ctselect/Makefile])
