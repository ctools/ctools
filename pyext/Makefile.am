## Process this file with automake to produce Makefile.in
#
# The actual implemented logic is to use swig if swig is
# available on the system

if ENABLE_PYTHON_BINDING
if HAVE_PYTHON
if WITH_PYTHON_BINDING

# Set SWIG interface filename
SWIG_INTERFACE = ctatools.i

# Define the targets that should be built before anything else.
# These are honored by 'make all', 'make check' and 'make install'.
if HAVE_SWIG
  BUILT_SOURCES = ctatools_wrap.cpp ctatools.py
endif

# Set files to be removed in clean
if HAVE_SWIG
  CLEANFILES = ctatools_wrap.cpp ctatools.py
endif

all-local: ctatools_wrap.cpp ctatools.py
	$(PYTHON_BUILD_PREFIX) $(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
if HAVE_SWIG
ctatools_wrap.cpp ctatools.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -I$(prefix)/share/gammalib/swig -I$(GAMMALIB)/share/gammalib/swig -o ctatools_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
endif # HAVE_SWIG

endif # WITH_PYTHON_BINDING
endif # HAVE_PYTHON
endif # ENABLE_PYTHON_BINDING
