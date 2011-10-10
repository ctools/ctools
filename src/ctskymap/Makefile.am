## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = ctskymap

# Sources
sources = ctskymap.cpp

# Libtool library
noinst_LTLIBRARIES = libctskymap.la
libctskymap_la_SOURCES = $(sources)

# Program definition
ctskymap_SOURCES = main.cpp
ctskymap_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctskymap_CPPFLAGS = @CPPFLAGS@
ctskymap_LDADD = libctskymap.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctskymap.par
