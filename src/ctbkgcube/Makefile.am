## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctbkgcube

# Sources
sources = ctbkgcube.cpp

# Libtool library
noinst_LTLIBRARIES = libctbkgcube.la
libctbkgcube_la_SOURCES = $(sources)

# Program definition
ctbkgcube_SOURCES = main.cpp
ctbkgcube_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctbkgcube_CPPFLAGS = @CPPFLAGS@
ctbkgcube_LDADD = libctbkgcube.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctbkgcube.par
