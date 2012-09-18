#############################################################################
# Makefile for ctools Python extension module                               #
# ------------------------------------------------------------------------- #
# Copyright (C) 2011-2012 Juergen Knoedlseder                               #
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

if WITH_PYTHON

# Set SWIG interface filename
SWIG_INTERFACE = ctools.i

# Define the targets that should be built before anything else.
# These are honored by 'make all', 'make check' and 'make install'.
if HAVE_SWIG
  BUILT_SOURCES = ctools_wrap.cpp ctools.py
endif

# Set files to be removed in clean
if HAVE_SWIG
  CLEANFILES = ctools_wrap.cpp ctools.py
endif

all-local: build

# We build the python extension twice since the install step needs the result
# in a specific directory and the test make check has no access to the name
# of this specific directory. So we build one instance for checking, the other
# for installation. If we could get the name built by distutils we could avoid
# this.
#$(PYTHON) setup.py build --build-platlib build
build: ctools_wrap.cpp ctools.py
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build _ctools.so

# Rule for building Python wrappers using SWIG.
if HAVE_SWIG
ctools_wrap.cpp ctools.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -I$(prefix)/share/gammalib/swig -I$(GAMMALIB)/share/gammalib/swig -o ctools_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
endif # HAVE_SWIG

endif # WITH_PYTHON
