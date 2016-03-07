## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctedispcube

# Sources
sources = ctedispcube.cpp

# Libtool library
noinst_LTLIBRARIES = libctedispcube.la
libctedispcube_la_SOURCES = $(sources)
libctedispcube_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctedispcube_SOURCES = main.cpp
ctedispcube_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctedispcube_CPPFLAGS = @CPPFLAGS@
ctedispcube_LDADD = libctedispcube.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctedispcube.par
