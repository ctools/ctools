## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctexpcube

# Sources
sources = ctexpcube.cpp

# Libtool library
noinst_LTLIBRARIES = libctexpcube.la
libctexpcube_la_SOURCES = $(sources)
libctexpcube_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctexpcube_SOURCES = main.cpp
ctexpcube_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctexpcube_CPPFLAGS = @CPPFLAGS@
ctexpcube_LDADD = libctexpcube.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctexpcube.par
