## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctcubemask

# Sources
sources = ctcubemask.cpp

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctcubemask.la
libctcubemask_la_SOURCES = $(sources)
libctcubemask_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctcubemask_SOURCES = main.cpp
ctcubemask_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctcubemask_CPPFLAGS = @CPPFLAGS@
ctcubemask_LDADD = libctcubemask.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctcubemask.par
