#############################################################################
# Makefile.am for ctools unit tests                                         #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2017 Juergen Knoedlseder                               #
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
#############################################################################

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Optional python bindings
if WITH_PYTHON
  TEST_PYTHON_CTOOLS   = test_python_ctools.sh
  TEST_PYTHON_CSCRIPTS = test_python_cscripts.sh
  TEST_PYTHON_EXAMPLES = test_examples.py
  TEST_PYTHON_ENV      = PYTHONPATH=$(top_builddir):$(top_builddir)/pyext:$(top_builddir)/pyext/build/ctools/ctools:$(pythondir):$(pyexecdir):$(PYTHONPATH)
endif

# Set test environment
TESTS_ENVIRONMENT = @RUNSHARED@=$(top_builddir)/src/.libs$(TEST_ENV_DIR):$(LIBDIR_GAMMALIB):$(@RUNSHARED@) \
                    $(TEST_PYTHON_ENV) \
                    TEST_SRCDIR=$(abs_top_srcdir) \
                    TEST_BUILDDIR=$(abs_top_builddir) \
                    TEST_DATA=$(abs_top_srcdir)/test/data

# Tests (these are the tests to be executed)
TESTS = $(TEST_PYTHON_CTOOLS) \
        $(TEST_PYTHON_CSCRIPTS) \
        $(TEST_PYTHON_EXAMPLES)

# Create report directory
all-local:
	mkdir -p reports

# Files to be removed in build directory
clean-local:
	rm -rf $(top_builddir)/test/*.pyc
	rm -rf $(top_builddir)/test/__pycache__
	rm -rf $(top_builddir)/test/*.fits
	rm -rf $(top_builddir)/test/*.log
	rm -rf $(top_builddir)/test/*.xml
	rm -rf $(top_builddir)/test/*.eps
	rm -rf $(top_builddir)/test/*.dat
	rm -rf $(top_builddir)/test/*.reg
	rm -rf $(top_builddir)/test/reports/*.xml
	rm -rf $(top_builddir)/test/caldb
	rm -rf $(top_builddir)/test/pfiles
	rm -rf $(top_builddir)/test/prod3
	rm -rf $(top_builddir)/test/csroot2caldb_*
	rm -rf $(top_builddir)/test/csiactcopy_*
	rm -rf $(top_builddir)/test/.DS_Store

# Tests to be done on installed version
if WITH_PYTHON
installcheck-local:
	@export CTOOLS=$(DESTDIR)$(prefix); \
	. $(DESTDIR)$(prefix)/bin/ctools-init.sh; \
	python -c 'import ctools; ctools.test()'; \
	python -c 'import cscripts; cscripts.test()'
endif

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
