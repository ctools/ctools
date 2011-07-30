## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = ctlike

# Sources
sources = ctlike.cpp

# Libtool library
noinst_LTLIBRARIES = libctlike.la
libctlike_la_SOURCES = $(sources)

# Program definition
ctlike_SOURCES = main.cpp
ctlike_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctlike_CPPFLAGS = @CPPFLAGS@
ctlike_LDADD = libctlike.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctlike.par
