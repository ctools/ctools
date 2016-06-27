## Process this file with automake to produce Makefile.in

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Optional python bindings
if WITH_PYTHON
  TEST_CSCRIPTS        = test_cscripts.sh
  TEST_PYTHON_CTOOLS   = test_python_ctools.sh
  TEST_PYTHON_CSCRIPTS = test_python_cscripts.sh
  TEST_PYTHON_EXAMPLES = test_examples.sh
  TEST_PYTHON_ENV      = PYTHONPATH=$(top_builddir):$(top_builddir)/pyext:$(top_builddir)/pyext/build/ctools/ctools:$(pythondir):$(pyexecdir):$(PYTHONPATH)
endif

# Set test environment
TESTS_ENVIRONMENT = @RUNSHARED@=$(top_builddir)/src/.libs$(TEST_ENV_DIR):$(LIBDIR_GAMMALIB):$(@RUNSHARED@) \
                    $(TEST_PYTHON_ENV)

# Tests (these are the tests to be executed)
TESTS = $(TEST_PYTHON_CTOOLS) \
        $(TEST_PYTHON_CSCRIPTS) \
        $(TEST_PYTHON_EXAMPLES)

# Add Valgrind rule
valgrind:
	@if type valgrind >/dev/null 2>&1; then \
		for prog in $(check_PROGRAMS) ; do \
			valgrind --leak-check=full \
			--xml=yes --xml-file=tmp.xml ./$$prog; \
			sed '/<\/valgrindoutput>/d' tmp.xml > valgrind_$$prog.xml; \
			echo '</valgrindoutput>' >> valgrind_$$prog.xml; \
		done ; \
	else \
		echo "valgrind  not found in PATH."; \
	fi
