## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = ctobssim

# Sources
sources = ctobssim.cpp

# Libtool library
noinst_LTLIBRARIES = libctobssim.la
libctobssim_la_SOURCES = $(sources)

# Program definition
ctobssim_SOURCES = main.cpp
ctobssim_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctobssim_CPPFLAGS = @CPPFLAGS@
ctobssim_LDADD = libctobssim.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctobssim.par
