#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for cscripts used from Python.
#
# Copyright (C) 2016 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import gammalib
import cscripts
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import test_cscaldb


# ================== #
# Perform unit tests #
# ================== #
def test(installed=False):
    """
    Perform unit testing for Python interface.
    """
    # Allocate test suite container
    suites = gammalib.GTestSuites("cscripts Python module unit testing")

    # Allocate test suites and append them to the container
    suite_cscaldb = test_cscaldb.Test()

    # Setup unit tests
    suite_cscaldb.set()

    # Append tests to container
    suites.append(suite_cscaldb)

    # If we have an installed version then create a temporary directory and
    # copy over all information that is needed
    if installed:

        # Create temporary working directory
        import tempfile
        path = tempfile.mkdtemp()
        os.chdir(path)

        # Get test directory
        import inspect
        testdir = inspect.getfile(cscripts.tests)
        dirname = os.path.dirname(testdir)

        # Copy over test data and irf
        #os.system("cp -r %s %s" % (dirname+"/data", "data"))
        #os.system("cp -r %s %s" % (dirname+"/irf",  "irf"))

    # Set PFILES environment variable
    try:
        os.mkdir("pfiles")
    except:
        pass
    os.environ['PFILES'] = "pfiles"

    # Copy over pfiles
    if not installed:
        os.system("cp -r ../cscripts/*.par pfiles/")
    else:
        os.system("cp -r %s/syspfiles/*.par pfiles/" % (os.environ['CTOOLS']))

    # Run test suite
    success = suites.run()

    # If we have a non-installed version then save test results
    if not installed:
        suites.save("reports/cscripts.xml")

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Exit with return code
    sys.exit(rc)


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Perform unit testing for cscripts used from Python.
    """
    # Run tests
    test()
