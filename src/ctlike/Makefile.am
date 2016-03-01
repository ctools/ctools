## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctlike

# Sources
sources = ctlike.cpp

# Libtool library
noinst_LTLIBRARIES = libctlike.la
libctlike_la_SOURCES = $(sources)
libctlike_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctlike_SOURCES = main.cpp
ctlike_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctlike_CPPFLAGS = @CPPFLAGS@
ctlike_LDADD = libctlike.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctlike.par
