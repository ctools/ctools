## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src

bin_PROGRAMS = ctbin

ctbin_SOURCES = ctbin.cpp
ctbin_LDFLAGS = @LDFLAGS@
ctbin_CPPFLAGS = @CPPFLAGS@

pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctbin.par
