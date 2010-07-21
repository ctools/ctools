## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

# Include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src

# Program name
bin_PROGRAMS = ctlike

# Sources
sources = ctlike.cpp

# Program definition
ctlike_SOURCES = main.cpp $(sources)
ctlike_LDFLAGS = @LDFLAGS@
ctlike_CPPFLAGS = @CPPFLAGS@

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctlike.par

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctlike.la
libctlike_la_SOURCES = $(sources)
