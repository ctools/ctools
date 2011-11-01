## Process this file with automake to produce Makefile.in

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Optional python bindings
if WITH_PYTHON
  TEST_CSCRIPTS      = test_cscripts.sh
  TEST_PYTHON_SCRIPT = test_python.py
  TEST_PYTHON_ENV    = PYTHONPATH=$(top_builddir)/pyext:$(top_builddir)/scripts:$(pythondir):$(pyexecdir):$(PYTHONPATH)
endif

# Set test environment
TESTS_ENVIRONMENT = $(TEST_PYTHON_ENV)

# Tests (these are the tests to be executed)
TESTS = test_ctools.sh \
        $(TEST_CSCRIPTS) \
        $(TEST_PYTHON_SCRIPT)
