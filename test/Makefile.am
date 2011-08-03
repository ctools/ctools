## Process this file with automake to produce Makefile.in

# Optional python bindings
if WITH_PYTHON
  TEST_CSCRIPTS      = test_cscripts.sh
  TEST_PYTHON_SCRIPT = test_python.py
  TEST_PYTHON_ENV    = PYTHONPATH=$(top_builddir)/pyext/build:$(top_builddir)/scripts:$(pythondir):$(PYTHONPATH)
endif

# Set test environment
TESTS_ENVIRONMENT = LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(prefix)/lib:$(LD_LIBRARY_PATH) \
                    PATH=$(top_builddir)/src/ctobssim:$(top_builddir)/src/ctbin:$(top_builddir)/src/ctselect:$(top_builddir)/src/ctlike:$(PATH) \
                    $(TEST_PYTHON_ENV)

# Tests (these are the tests to be executed)
TESTS = test_ctools.sh \
        $(TEST_CSCRIPTS) \
        $(TEST_PYTHON_SCRIPT)
