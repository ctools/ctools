## Process this file with automake to produce Makefile.in

# Set sub directories
SUBDIRS = ctselect \
          ctbin \
          ctlike \
          ctobssim \
          ctskymap

# Scripts to be distributed
scripts = ctools-setup \
          ctools-init.sh \
          ctools-init.csh

# Build libtool library for ctools. This is needed to provide
# a shared library that can be loaded for python bindings. We do
# not have to build this library if no python bindings are used.
lib_LTLIBRARIES = libctools.la
libctools_la_SOURCES = 
libctools_la_LIBADD = ctbin/libctbin.la \
                      ctselect/libctselect.la \
                      ctlike/libctlike.la \
                      ctobssim/libctobssim.la \
                      ctskymap/libctskymap.la

# Distribute scripts
EXTRA_DIST   = $(scripts)
scriptsdir   = $(bindir)
scripts_DATA = $(scripts)

# Add hook rule to make scripts executable
install-data-hook:
	@list='$(scripts_DATA)';\
	for p in $$list; do \
	  if test -f "$$p"; then chmod +x "$(DESTDIR)$(scriptsdir)/$$p"; fi; \
	done
