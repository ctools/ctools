## Process this file with automake to produce Makefile.in
#
# The actual implemented logic is to use swig if swig is
# available on the system

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
build: ctools_wrap.cpp ctools.py
	$(PYTHON_BUILD_PREFIX) $(PYTHON) setup.py build --build-platlib build
	$(PYTHON_BUILD_PREFIX) $(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
if HAVE_SWIG
ctools_wrap.cpp ctools.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -I$(prefix)/share/gammalib/swig -I$(GAMMALIB)/share/gammalib/swig -o ctools_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
endif # HAVE_SWIG

endif # WITH_PYTHON
