#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for cscripts
#
# Copyright (C) 2016-2017 Juergen Knoedlseder
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
import test_cscript
import test_cscaldb
import test_csebins
import test_csfindobs
import test_csinfo
import test_cslightcrv
import test_csmodelinfo
import test_csmodelmerge
import test_csmodelselect
import test_csobs2caldb
import test_csobsdef
import test_csobsinfo
import test_csobsselect
import test_cspull
import test_csresmap
import test_csroot2caldb
import test_cssens
import test_csspec
import test_cssrcdetect
import test_cstsdist
import test_cstsmapmerge
import test_cstsmapsplit
import test_csviscube
import test_csworkflow
import test_csiactcopy
import test_csiactdata
import test_csiactobs


# ================== #
# Perform unit tests #
# ================== #
def test(installed=False, debug=False):
    """
    Perform cscripts unit testing.

    Parameters
    ----------
    installed : bool, optional
        Flag indicating whether the script has been installed or not
    debug : bool, optional
        Enables debugging
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

        # Copy test data in "data" directory
        os.system('cp -r %s %s' % (dirname+'/data', 'data'))

        # Set test data environment variable
        os.environ['TEST_DATA'] = 'data'

    # ... otherwise set the calibration database to the one shipped with the
    # package; we don't need to set the 'TEST_DATA', this is done by the
    # test environment
    else:
        os.environ['CALDB'] = '%s/caldb' % (os.environ['TEST_SRCDIR'])

    # Create a local "pfiles" directory and set PFILES environment variable
    try:
        os.mkdir('pfiles')
    except:
        pass
    os.environ['PFILES'] = 'pfiles'

    # Copy the cscripts parameter files into the "pfiles" directory. For a
    # non-installed test we copy the parameter files from the "cscripts"
    # source directory into the "pfiles" directory, for an installed version
    # we get all parameter files from the "syspfiles" directory. Also make
    # sure that all parameter files are writable.
    if not installed:
        os.system('cp -r %s/cscripts/*.par pfiles/'% (os.environ['TEST_SRCDIR']))
        os.system('chmod u+w pfiles/*')
    else:
        os.system('cp -r %s/syspfiles/*.par pfiles/' % (os.environ['CTOOLS']))
        os.system('chmod u+w pfiles/*')

    # Define list of test suites
    tests = [test_cscript.Test(),
             test_cscaldb.Test(),
             test_csebins.Test(),
             test_csinfo.Test(),
             test_cslightcrv.Test(),
             test_csmodelinfo.Test(),
             test_csmodelmerge.Test(),
             test_csmodelselect.Test(),
             test_csobs2caldb.Test(),
             test_csobsdef.Test(),
             test_csobsinfo.Test(),
             test_csobsselect.Test(),
             test_cspull.Test(),
             test_csresmap.Test(),
             test_csroot2caldb.Test(),
             test_cssens.Test(),
             test_csspec.Test(),
             test_cssrcdetect.Test(),
             test_cstsdist.Test(),
             test_cstsmapmerge.Test(),
             test_cstsmapsplit.Test(),
             test_csviscube.Test(),
             test_csworkflow.Test()]

    # Append tests for Python 2.6+ (the IACT cscripts depend on the json
    # module which is only available since Python 2.6+)
    ver = sys.version.split()[0]
    if ver >= '2.6.0':

        # Check for VHEFITS environment variable
        if 'VHEFITS' in os.environ:
            
            # If the environment variable exists then unset it for the test
            # cases. Since Python is executed in a subprocess this will not
            # impact the environment variable in the parent shell.
            del os.environ['VHEFITS']

        # Append IACT script tests to list of test suites
        tests.extend([test_csfindobs.Test(),
                      test_csiactcopy.Test(),
                      test_csiactdata.Test(),
                      test_csiactobs.Test()])

    # Allocate test suite container
    suites = gammalib.GTestSuites('cscripts unit testing')

    # Set test suites and append them to suite container
    for suite in tests:
        suite.set()
        suites.append(suite)

    # Run test suite
    success = suites.run()

    # If we have a non-installed version then save test results
    if not installed:
        suites.save('reports/cscripts.xml')

    # If debuging is requested then print test suites
    if debug:
        print(suites)

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Remove temporary direction
    if installed:
        os.system('rm -rf %s' % (path))

    # Exit with return code
    sys.exit(rc)


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run tests
    test()
