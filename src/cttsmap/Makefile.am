## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = cttsmap

# Sources
sources = cttsmap.cpp

# Libtool library
noinst_LTLIBRARIES = libcttsmap.la
libcttsmap_la_SOURCES = $(sources)
libcttsmap_la_LDFLAGS = -version-info $(CTOOLS_LT_VERSION) @OPENMP_LIB_LDFLAGS@

# Program definition
cttsmap_SOURCES = main.cpp
cttsmap_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
cttsmap_CPPFLAGS = @CPPFLAGS@
cttsmap_LDADD = libcttsmap.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = cttsmap.par
