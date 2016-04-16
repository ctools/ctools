## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctmapcube

# Sources
sources = ctmapcube.cpp

# Libtool library
noinst_LTLIBRARIES = libctmapcube.la
libctmapcube_la_SOURCES = $(sources)
libctmapcube_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctmapcube_SOURCES = main.cpp
ctmapcube_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctmapcube_CPPFLAGS = @CPPFLAGS@
ctmapcube_LDADD = libctmapcube.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctmapcube.par
