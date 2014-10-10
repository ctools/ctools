## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctbin

# Sources
sources = ctbin.cpp

# Libtool library
noinst_LTLIBRARIES = libctbin.la
libctbin_la_SOURCES = $(sources)

# Program definition
ctbin_SOURCES = main.cpp
ctbin_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctbin_CPPFLAGS = @CPPFLAGS@
ctbin_LDADD = libctbin.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctbin.par
