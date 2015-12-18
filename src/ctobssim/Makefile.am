## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctobssim

# Sources
sources = ctobssim.cpp

# Libtool library
noinst_LTLIBRARIES = libctobssim.la
libctobssim_la_SOURCES = $(sources)
libctobssim_la_LDFLAGS = -version-info $(CTOOLS_LT_VERSION) @OPENMP_LIB_LDFLAGS@

# Program definition
ctobssim_SOURCES = main.cpp
ctobssim_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctobssim_CPPFLAGS = @CPPFLAGS@
ctobssim_LDADD = libctobssim.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctobssim.par
