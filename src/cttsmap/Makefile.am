## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = cttsmap

# Sources
sources = cttsmap.cpp

# Libtool library
noinst_LTLIBRARIES = libcttsmap.la
libcttsmap_la_SOURCES = $(sources)

# Program definition
cttsmap_SOURCES = main.cpp
cttsmap_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
cttsmap_CPPFLAGS = @CPPFLAGS@
cttsmap_LDADD = libcttsmap.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = cttsmap.par
