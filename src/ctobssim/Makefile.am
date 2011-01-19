## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

# Include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src

# Program name
bin_PROGRAMS = ctobssim

# Sources
sources = ctobssim.cpp

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctobssim.la
libctobssim_la_SOURCES = $(sources)

# Program definition
ctobssim_SOURCES = main.cpp $(sources)
ctobssim_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctobssim_CPPFLAGS = @CPPFLAGS@

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctobssim.par
