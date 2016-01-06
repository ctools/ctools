## Process this file with automake to produce Makefile.in

# Define sources for this directory
sources = ctool.cpp

# Build libtool library
noinst_LTLIBRARIES = libsupport.la
libsupport_la_SOURCES = $(sources)
libsupport_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@
