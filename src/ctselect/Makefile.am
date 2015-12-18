## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctselect

# Sources
sources = ctselect.cpp

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctselect.la
libctselect_la_SOURCES = $(sources)
libctselect_la_LDFLAGS = -version-info $(CTOOLS_LT_VERSION) @OPENMP_LIB_LDFLAGS@

# Program definition
ctselect_SOURCES = main.cpp
ctselect_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctselect_CPPFLAGS = @CPPFLAGS@
ctselect_LDADD = libctselect.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctselect.par
