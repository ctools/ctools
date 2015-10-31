## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = cterror

# Sources
sources = cterror.cpp

# Libtool library
noinst_LTLIBRARIES = libcterror.la
libcterror_la_SOURCES = $(sources)

# Program definition
cterror_SOURCES = main.cpp
cterror_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
cterror_CPPFLAGS = @CPPFLAGS@
cterror_LDADD = libcterror.la $(top_srcdir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = cterror.par
