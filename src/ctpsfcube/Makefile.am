## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctpsfcube

# Sources
sources = ctpsfcube.cpp

# Libtool library
noinst_LTLIBRARIES = libctpsfcube.la
libctpsfcube_la_SOURCES = $(sources)
libctpsfcube_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctpsfcube_SOURCES = main.cpp
ctpsfcube_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctpsfcube_CPPFLAGS = @CPPFLAGS@
ctpsfcube_LDADD = libctpsfcube.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctpsfcube.par
