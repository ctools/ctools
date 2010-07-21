## Process this file with automake to produce Makefile.in

if HAVE_SWIG
if HAVE_PYTHON
if ENABLE_PYTHON_BINDING

SWIG_INTERFACE = ctatools.i
BUILT_SOURCES = ctatools_wrap.cpp ctatools.py
CLEANFILES = ctatools_wrap.cpp ctatools.py

all-local: ctatools_wrap.cpp ctatools.py
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
ctatools_wrap.cpp ctatools.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -includeall -I$(prefix)/share/gammalib/swig -o ctatools_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
	@if test "x`uname`" = "xCygwin"; then \
		echo "WARNING: working around bug in SWIG on Windows/Cygwin platforms ONLY"; \
		sed -i.bak -f $(srcdir)/fix_wrapper.sed ctatools_wrap.cpp; \
	fi

endif # ENABLE_PYTHON_BINDING
endif # HAVE_PYTHON
endif # HAVE_SWIG
