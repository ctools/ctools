## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctskymap

# Sources
sources = ctskymap.cpp

# Libtool library
noinst_LTLIBRARIES = libctskymap.la
libctskymap_la_SOURCES = $(sources)
libctskymap_la_LDFLAGS = -version-info $(CTOOLS_LT_VERSION) @OPENMP_LIB_LDFLAGS@

# Program definition
ctskymap_SOURCES = main.cpp
ctskymap_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctskymap_CPPFLAGS = @CPPFLAGS@
ctskymap_LDADD = libctskymap.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctskymap.par
