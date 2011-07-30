## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = ctselect

# Sources
sources = ctselect.cpp

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctselect.la
libctselect_la_SOURCES = $(sources)

# Program definition
ctselect_SOURCES = main.cpp
ctselect_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctselect_CPPFLAGS = @CPPFLAGS@
ctselect_LDADD = libctselect.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctselect.par
