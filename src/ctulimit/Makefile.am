## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctulimit

# Sources
sources = ctulimit.cpp

# Libtool library
noinst_LTLIBRARIES = libctulimit.la
libctulimit_la_SOURCES = $(sources)
libctulimit_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctulimit_SOURCES = main.cpp
ctulimit_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctulimit_CPPFLAGS = @CPPFLAGS@
ctulimit_LDADD = libctulimit.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctulimit.par
