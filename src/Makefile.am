## Process this file with automake to produce Makefile.in

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Set sub directories
SUBDIRS = support \
          ctbin \
          ctlike \
          ctmodel \
          ctobssim \
          ctselect \
          ctskymap \
          ctexpcube \
          ctpsfcube \
          ctedispcube \
          ctbkgcube \
          ctcubemask \
          cttsmap \
          ctbutterfly \
          ctulimit \
          cterror

# Scripts to be distributed
scripts = ctools-setup \
          ctools-init.sh \
          ctools-init.csh

# Build libtool library for ctools. This is needed to provide
# a shared library that can be loaded for python bindings. We do
# not have to build this library if no python bindings are used.
lib_LTLIBRARIES = libctools.la
libctools_la_SOURCES = 
libctools_la_LIBADD = support/libsupport.la \
                      ctbin/libctbin.la \
                      ctlike/libctlike.la \
                      ctmodel/libctmodel.la \
                      ctobssim/libctobssim.la \
                      ctselect/libctselect.la \
                      ctskymap/libctskymap.la \
                      ctexpcube/libctexpcube.la \
                      ctpsfcube/libctpsfcube.la \
                      ctedispcube/libctedispcube.la \
                      ctbkgcube/libctbkgcube.la \
                      ctcubemask/libctcubemask.la \
                      cttsmap/libcttsmap.la \
                      ctbutterfly/libctbutterfly.la \
                      ctulimit/libctulimit.la \
                      cterror/libcterror.la
libctools_la_LDFLAGS = -version-info $(CTOOLS_LT_VERSION) @OPENMP_LIB_LDFLAGS@

# Distribute scripts
EXTRA_DIST   = $(scripts)
scriptsdir   = $(bindir)
scripts_DATA = $(scripts)

# Things to install
install-data-local: help-install

# Install help files
help-install:
	@$(mkinstalldirs) $(DESTDIR)$(datadir)/help && \
    for file in $(top_srcdir)/doc/source/reference_manual/c*.rst ; do \
        filename="$${file##*/}" ; \
        basename="$${filename%%.*}" ; \
        helpname="$$basename.txt" ; \
        echo "$(INSTALL_DATA) $$file $(DESTDIR)$(datadir)/help/$$helpname" ; \
        $(INSTALL_DATA) "$$file" $(DESTDIR)$(datadir)/help/"$$helpname" ; \
        sed -i.bak -e 's/``//g' -e 's/`//g' -e 's/:ref://g' -e 's/:doc://g' \
            -e 's/.. note::/Note/g' -e '1d' \
            $(DESTDIR)$(datadir)/help/"$$helpname" ; \
        rm -rf $(DESTDIR)$(datadir)/help/"$$helpname".bak ; \
    done

# Add hook rule to make scripts executable
install-data-hook:
	@list='$(scripts_DATA)';\
	for p in $$list; do \
	  if test -f "$$p"; then chmod +x "$(DESTDIR)$(scriptsdir)/$$p"; fi; \
	done
