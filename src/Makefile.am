## Process this file with automake to produce Makefile.in

# Set sub directories
SUBDIRS = ctselect ctbin ctlike

# Build libtool library for CTA tools. This is needed to provide
# a shared library that can be loaded for python bindings. We do
# not have to build this library if no python bindings are used.
lib_LTLIBRARIES = libctatools.la
libctatools_la_SOURCES = 
libctatools_la_LIBADD = ctbin/libctbin.la \
		        ctselect/libctselect.la \
                        ctlike/libctlike.la
