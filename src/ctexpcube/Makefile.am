## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = ctexpcube

# Sources
sources = ctexpcube.cpp

# Libtool library
noinst_LTLIBRARIES = libctexpcube.la
libctexpcube_la_SOURCES = $(sources)

# Program definition
ctexpcube_SOURCES = main.cpp
ctexpcube_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctexpcube_CPPFLAGS = @CPPFLAGS@
ctexpcube_LDADD = libctexpcube.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctexpcube.par
