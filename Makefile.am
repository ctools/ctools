## Process this file with automake to produce Makefile.in

# Include M4 macros
ACLOCAL_AMFLAGS = -I m4

# Subdirectories
SUBDIRS = . src pyext test

#
INCLUDES =

# Make sure these will be cleaned even when they're not built by default
CLEANFILES =

# Headers to be installed
include_HEADERS =

# Headers which are not installed but which are distributed
noinst_HEADERS =

# Pkgconfig directory
pkgconfigdir = $(libdir)/pkgconfig

# Files to install in Pkgconfig directory
pkgconfig_DATA = ctatools.pc

# Manual pages to install
man_MANS = 

# Architecture independent data files installed in the package's doc directory
doc_DATA = \
	   AUTHORS \
	   ChangeLog \
	   COPYING \
	   INSTALL \
	   LICENSE \
	   NEWS \
	   README
docdir = ${datadir}/doc/${PACKAGE}

# Models installed in the package's model directory
models      = test/data/crab.xml
modelsdir   = $(datadir)/models
models_DATA = $(models)

# Models to install
# Extra files to be distributed in the tar file
#EXTRA_DIST = $(models)


install-data-local: $(HTML_INSTALL_DATA_TARGETS)

uninstall-local: $(HTML_UNINSTALL_DATA_TARGETS)


# Non-Automake subdirectories to distribute
DISTDIRS = doc/html doc/html/images
dist-hook:
	( \
	  builddir=`pwd` ; \
	  cd $(srcdir) && \
	  ( \
	    for dir in $(DISTDIRS) ; do \
	      find $$dir -depth -print | egrep -v '(~$$)|(/CVS)|(/\.#)|(/\.deps)' \
		| cpio -pdum $$builddir/$(distdir) 2> /dev/null ; \
	    done \
	  ) \
	)


#
# Additional install rules
#

# Install HTML files
DOCDIRS = html html/images
HTML_INSTALL_DATA_TARGETS = install-data-html
install-data-html:
	$(mkinstalldirs) $(DESTDIR)$(docdir)
	@for dir in $(DOCDIRS) ; do \
	  $(mkinstalldirs) $(DESTDIR)$(docdir)/$$dir && \
	  for file in $(srcdir)/doc/$$dir/*.* ; do \
	    echo "$(INSTALL_DATA) $$file $(DESTDIR)$(docdir)/$$dir" ; \
	    $(INSTALL_DATA) "$$file" $(DESTDIR)$(docdir)/$$dir ; \
	  done ; \
	done

# Uninstall HTML files
HTML_UNINSTALL_DATA_TARGETS = uninstall-data-html
uninstall-data-html:
	for dir in $(DOCDIRS) ; do \
	  rm -f -r $(DESTDIR)$(docdir)/$$dir ; \
	done

# Automatically reconfigure libtool
LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck
