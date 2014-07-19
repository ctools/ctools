## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = ctcubemask

# Sources
sources = ctcubemask.cpp

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctcubemask.la
libctcubemask_la_SOURCES = $(sources)

# Program definition
ctcubemask_SOURCES = main.cpp
ctcubemask_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctcubemask_CPPFLAGS = @CPPFLAGS@
ctcubemask_LDADD = libctcubemask.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctcubemask.par
