## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

# Include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src

# Program name
bin_PROGRAMS = ctbin

# Sources
sources = ctbin.cpp

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctbin.la
libctbin_la_SOURCES = $(sources)

# Program definition
ctbin_SOURCES = main.cpp $(sources)
ctbin_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctbin_CPPFLAGS = @CPPFLAGS@

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctbin.par
