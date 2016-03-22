#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for cscripts.
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
import os
import sys
import gammalib
import cscripts
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import test_cscaldb
import test_csobs2caldb
import test_cslightcrv
import test_csfindobs
import test_csiactcopy
import test_csiactdata


# ================== #
# Perform unit tests #
# ================== #
def test(installed=False, debug=False):
    """
    Perform cscripts unit testing.
    """
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
        os.system("cp -r %s %s" % (dirname+"/data", "data"))
        os.system("cp -r %s %s" % (dirname+"/irf",  "irf"))

    # ... otherwise set the calibration database
    else:
        os.environ['CALDB'] = "../caldb"

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

    # Allocate test suite container
    suites = gammalib.GTestSuites("cscripts unit testing")

    # Allocate test suites
    suite_cscaldb     = test_cscaldb.Test()
    suite_csobs2caldb = test_csobs2caldb.Test()
    suite_cslightcrv  = test_cslightcrv.Test()

    # Setup unit tests
    suite_cscaldb.set()
    suite_csobs2caldb.set()
    suite_cslightcrv.set()

    # Append tests to container
    suites.append(suite_cscaldb)
    suites.append(suite_csobs2caldb)
    suites.append(suite_cslightcrv)

    # Append tests for Python 2.6+ (the IACT cscripts depend on the json
    # module which is only available since Python 2.6+
    ver = sys.version.split()[0]
    if ver >= '2.6.0':

        # Allocate test suites
        suite_csfindobs  = test_csfindobs.Test()
        suite_csiactcopy = test_csiactcopy.Test()
        suite_csiactdata = test_csiactdata.Test()

        # Setup unit tests
        suite_csfindobs.set()
        suite_csiactcopy.set()
        suite_csiactdata.set()

        # Append tests to container
        suites.append(suite_csfindobs)
        suites.append(suite_csiactcopy)
        suites.append(suite_csiactdata)

    # Run test suite
    success = suites.run()

    # If we have a non-installed version then save test results
    if not installed:
        suites.save("reports/cscripts.xml")

    # If debuging is requested then print test suites
    if debug:
        print(suites)

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

    # Run tests
    test()
