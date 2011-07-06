## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

# Optional python bindings
if WITH_PYTHON
  TEST_PYTHON_SCRIPT = test_python.py
  TEST_PYTHON_ENV    = PYTHONPATH=$(top_builddir)/pyext/build:$(pythondir):$(PYTHONPATH)
endif

# Include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src

# Test programs (those will be compiled)
check_PROGRAMS = 

# Set test environment
TESTS_ENVIRONMENT = LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(prefix)/lib:$(LD_LIBRARY_PATH) \
                    PATH=$(top_builddir)/src/ctobssim:$(top_builddir)/src/ctbin:$(top_builddir)/src/ctselect:$(top_builddir)/src/ctlike:$(PATH) \
                    $(TEST_PYTHON_ENV)

# Tests (these are the tests to be executed)
TESTS = test_ctatools.sh \
        $(TEST_PYTHON_SCRIPT)

# Standard test sources and links
#test_GVector_SOURCES = test_GVector.cpp
#test_GVector_LDFLAGS = @LDFLAGS@
#test_GVector_CPPFLAGS = @CPPFLAGS@
#test_GVector_LDADD = $(top_srcdir)/src/libgamma.la
