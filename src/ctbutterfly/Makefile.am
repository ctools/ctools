## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctbutterfly

# Sources
sources = ctbutterfly.cpp

# Libtool library
noinst_LTLIBRARIES = libctbutterfly.la
libctbutterfly_la_SOURCES = $(sources)

# Program definition
ctbutterfly_SOURCES = main.cpp
ctbutterfly_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctbutterfly_CPPFLAGS = @CPPFLAGS@
ctbutterfly_LDADD = libctbutterfly.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctbutterfly.par
