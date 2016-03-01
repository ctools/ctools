## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctmodel

# Sources
sources = ctmodel.cpp

# Libtool library
noinst_LTLIBRARIES = libctmodel.la
libctmodel_la_SOURCES = $(sources)
libctmodel_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctmodel_SOURCES = main.cpp
ctmodel_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctmodel_CPPFLAGS = @CPPFLAGS@
ctmodel_LDADD = libctmodel.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctmodel.par
