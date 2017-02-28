#############################################################################
# Makefile for ctools Python extension module                               #
# ------------------------------------------------------------------------- #
# Copyright (C) 2011-2017 Juergen Knoedlseder                               #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
# ------------------------------------------------------------------------- #
# This Makefile will build the Python interface files (*.py) and the Python #
# wrapper files using SWIG (if available; otherwise we have to make sure    #
# the the files exist already). The Python distutils are used to build and  #
# install the Python modules (see file setup.py.in).                        #
#                                                                           #
# Process this file with automake to produce Makefile.in.                   #
#############################################################################

# Works only if Python is available
if WITH_PYTHON

# Set wrapper directory depending on the whether the wrappers were actually
# shipped with the distribution of whether the wrappers have to be build
# using swig
if BUILD_WRAPPERS
  wrapperdir = $(top_builddir)/pyext
else
  wrapperdir = $(srcdir)
endif

# Test scripts to be distributed in $(pyexecdir)/ctools/tests
testscripts = ctools/tests/__init__.py \
              $(top_srcdir)/test/testing.py \
              $(top_srcdir)/test/test_python_ctools.py \
              $(top_srcdir)/test/test_ctobssim.py \
              $(top_srcdir)/test/test_ctselect.py \
              $(top_srcdir)/test/test_ctbin.py \
              $(top_srcdir)/test/test_ctlike.py \
              $(top_srcdir)/test/test_cttsmap.py \
              $(top_srcdir)/test/test_ctmodel.py \
              $(top_srcdir)/test/test_ctskymap.py \
              $(top_srcdir)/test/test_ctexpcube.py \
              $(top_srcdir)/test/test_ctpsfcube.py \
              $(top_srcdir)/test/test_ctedispcube.py \
              $(top_srcdir)/test/test_ctbkgcube.py \
              $(top_srcdir)/test/test_ctmapcube.py \
              $(top_srcdir)/test/test_ctcubemask.py \
              $(top_srcdir)/test/test_ctbutterfly.py \
              $(top_srcdir)/test/test_ctulimit.py \
              $(top_srcdir)/test/test_cterror.py \
              $(top_srcdir)/test/test_pipelines.py

# Test data to be distributed in $(pyexecdir)/ctools/tests/data
testdata = $(top_srcdir)/test/data/crab.xml \
           $(top_srcdir)/test/data/crab_events.fits.gz \
           $(top_srcdir)/test/data/crab_cntmap.fits.gz \
           $(top_srcdir)/test/data/invalid_event_list.fits.gz \
           $(top_srcdir)/test/data/exclusion.reg

# Install test scripts and test data. Note that we don't need to prepend
# $(DESTDIR) here (not fully clear why, but that's the way it is)

testscriptsdir   = $(pyexecdir)/ctools/tests
testscripts_DATA = $(testscripts)
testdatadir      = $(pyexecdir)/ctools/tests/data
testdata_DATA    = $(testdata)

# If building of Python wrappers was requested then add them to the build
# source
if BUILD_WRAPPERS

# Set SWIG compile command
SWIGCOMPILE = $(SWIG) -c++ -python -Wall \
              -I$(DESTDIR)$(prefix)/share/gammalib/swig \
              -I$(GAMMALIB)/share/gammalib/swig \
              -I$(top_srcdir)/pyext

# Sources that are required for the build. All targets in the list will be
# built before anything. These are honored by 'make all', 'make check' and
# 'make install'. We use this here to make sure that all these targets are
# build by swig before we start compiling the Python modules.
BUILT_SOURCES = $(wrapperdir)/ctools/tools_wrap.cpp \
                $(wrapperdir)/ctools/tools.py

# Add the wrappers to the files that will be deleted in a 'make clean'
CLEAN_FILES = $(BUILT_SOURCES)

# Rule to build SWIG wrappers
swig: $(BUILT_SOURCES)

# Rule for building Python wrappers using SWIG.
$(wrapperdir)/ctools/tools.py: $(wrapperdir)/ctools/tools_wrap.cpp
$(wrapperdir)/ctools/tools_wrap.cpp: ctools/tools.i
	if $(SWIGCOMPILE) -MMD -MF "tools.Tpi" -o ctools/tools_wrap.cpp $<; \
	then mv -f "tools.Tpi" "tools.Pi"; else rm -f "tools.Tpi"; exit 1; fi

# ... otherwise, if building of Python wrappers has not been required then
# they are located in the source directory and we don't have to build them
# using swig
else
CLEAN_FILES =
BUILT_SOURCES =
endif # BUILD_WRAPPERS

# Extra files to be distributed in the tar file
EXTRA_DIST = ctools/tests \
             $(wrapperdir)/ctools/tools_wrap.cpp \
             $(wrapperdir)/ctools/tools.py \
             $(testscripts) \
             $(testdata) \
             $(testirf)

# Dependencies. Obtain the list of dependencies by executing ls.
DEPFILES = $(shell ls *.Pi 2> /dev/null)
include $(DEPFILES)

# Things we want to do when 'make' is invoked.
all-local: build

# Build the ctools extension module
build: $(BUILT_SOURCES)
	$(PYTHON) setup.py build
if IS_ELCAPITAN
	-for f in build/ctools/ctools/_*.so; do \
		echo "Build post process Python module "$$f; \
		libgamma=`otool -L $$f | grep -o libgamma.*.dylib`; \
		libctools=`otool -L $$f | grep -o libctools.*.dylib`; \
		install_name_tool -change $(libdir)/$$libctools @rpath/$$libctools $$f; \
		install_name_tool -change $(libdir)/$$libgamma @rpath/$$libgamma $$f; \
		install_name_tool -add_rpath $(abs_top_builddir)/src/.libs $$f; \
		install_name_tool -add_rpath ${GAMMALIB}/lib $$f; \
	done
endif

# Install the ctools Python module
install-exec-local:
if IS_ELCAPITAN
	-for f in build/ctools/ctools/_*.so; do \
		echo "Install pre process Python module "$$f; \
		install_name_tool -rpath $(abs_top_builddir)/src/.libs $(libdir) $$f; \
	done
endif
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

# Uninstall the ctools Python module
uninstall-local:
	rm -rf $(DESTDIR)$(pyexecdir)/ctools/*.py
	rm -rf $(DESTDIR)$(pyexecdir)/ctools/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/ctools/*.so
	rm -rf $(DESTDIR)$(pyexecdir)/ctools/__pycache__
	rm -rf $(DESTDIR)$(pyexecdir)/ctools/tests/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/ctools/tests/__pycache__
	rm -rf $(DESTDIR)$(pyexecdir)/ctools-@VERSION@*.egg-info

# Set files to be removed in 'make clean'
CLEANFILES = $(CLEAN_FILES) $(DEPFILES)

# Files to be removed in build directory
clean-local:
	rm -rf build
	rm -rf $(top_builddir)/pyext/ctools/tools.py
	rm -rf $(top_builddir)/pyext/ctools/c*.py
	rm -rf $(top_builddir)/pyext/ctools/*.pyc
	rm -rf $(top_builddir)/pyext/ctools/*_wrap.cpp
	rm -rf $(top_builddir)/pyext/ctools/_tools.so
	rm -rf $(top_builddir)/pyext/ctools/__pycache__

endif # WITH_PYTHON
