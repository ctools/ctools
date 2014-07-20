## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = ctpsfcube

# Sources
sources = ctpsfcube.cpp

# Libtool library
noinst_LTLIBRARIES = libctpsfcube.la
libctpsfcube_la_SOURCES = $(sources)

# Program definition
ctpsfcube_SOURCES = main.cpp
ctpsfcube_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctpsfcube_CPPFLAGS = @CPPFLAGS@
ctpsfcube_LDADD = libctpsfcube.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctpsfcube.par
