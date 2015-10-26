#############################################################################
# Makefile for ctools Python extension module                               #
# ------------------------------------------------------------------------- #
# Copyright (C) 2011-2015 Juergen Knoedlseder                               #
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

# Test scripts to be distributed in $(pyexecdir)/ctools/tests
testscripts = ctools/tests/__init__.py \
              $(abs_top_srcdir)/test/test_python.py \
              $(abs_top_srcdir)/test/test_ctobssim.py \
              $(abs_top_srcdir)/test/test_ctselect.py \
              $(abs_top_srcdir)/test/test_ctbin.py \
              $(abs_top_srcdir)/test/test_ctlike.py \
              $(abs_top_srcdir)/test/test_cttsmap.py \
              $(abs_top_srcdir)/test/test_ctmodel.py \
              $(abs_top_srcdir)/test/test_ctskymap.py \
              $(abs_top_srcdir)/test/test_ctexpcube.py \
              $(abs_top_srcdir)/test/test_ctpsfcube.py \
              $(abs_top_srcdir)/test/test_ctbkgcube.py \
              $(abs_top_srcdir)/test/test_ctcubemask.py \
              $(abs_top_srcdir)/test/test_ctbutterfly.py \
              $(abs_top_srcdir)/test/test_ctulimit.py \
              $(abs_top_srcdir)/test/test_cterror.py \
              $(abs_top_srcdir)/test/test_pipelines.py
testdata    = $(abs_top_srcdir)/test/data/crab.xml \
              $(abs_top_srcdir)/test/data/crab_events.fits.gz \
              $(abs_top_srcdir)/test/data/crab_cntmap.fits.gz \
              $(abs_top_srcdir)/test/data/exclusion.reg
testirf     = $(abs_top_srcdir)/test/irf/cta_dummy_irf.dat

# Distribute scripts
EXTRA_DIST       = $(testscripts) $(testdata) $(testirf)
testscriptsdir   = $(pyexecdir)/ctools/tests
testscripts_DATA = $(testscripts)
testdatadir      = $(pyexecdir)/ctools/tests/data
testdata_DATA    = $(testdata)
testirfdir       = $(pyexecdir)/ctools/tests/irf
testirf_DATA     = $(testirf)

# Set SWIG interface filename
SWIG_INTERFACE = ctools/tools.i

# Dependencies. Obtain the list of dependencies by executing ls.
DEPFILES = $(shell ls *.Pi 2> /dev/null)
include $(DEPFILES)

# Define the targets that should be built before anything else.
# These are honored by 'make all', 'make check' and 'make install'.
if HAVE_SWIG
  BUILT_SOURCES = ctools/tools_wrap.cpp ctools/tools.py
endif

# Set files to be removed in clean. We remove only the SWIG interface files
# if SWIG is available to rebuild them. 
if HAVE_SWIG
  CLEAN_FILES = $(BUILT_SOURCES)
endif
CLEANFILES = $(CLEAN_FILES) $(DEPFILES)

# Set SWIG compile command
SWIGCOMPILE = $(SWIG) -c++ -python -Wall -I$(prefix)/share/gammalib/swig -I$(GAMMALIB)/share/gammalib/swig

all-local: build

# We build the python extension twice since the install step needs the result
# in a specific directory and the test make check has no access to the name
# of this specific directory. So we build one instance for checking, the other
# for installation. If we could get the name built by distutils we could avoid
# this.
build: ctools/tools_wrap.cpp ctools/tools.py
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build _tools.so

# Rule for building Python wrappers using SWIG.
if HAVE_SWIG
ctools/tools.py: ctools/tools_wrap.cpp
ctools/tools_wrap.cpp: $(SWIG_INTERFACE)
	if $(SWIGCOMPILE) -MMD -MF "tools.Tpi" -o ctools/tools_wrap.cpp $<; \
	then mv -f "tools.Tpi" "tools.Pi"; else rm -f "tools.Tpi"; exit 1; fi
endif # HAVE_SWIG

endif # WITH_PYTHON
