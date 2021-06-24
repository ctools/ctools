#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for comscripts
#
# Copyright (C) 2021 Juergen Knoedlseder
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
import comscripts
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import test_comobsselect
import test_comobsbin
import test_comobsback
import test_comobsmodel
import test_comobsadd
import test_comobsres
import test_comobssim
import test_comlixfit
import test_comlixmap
import test_comsrcdetect


# ================== #
# Perform unit tests #
# ================== #
def test(installed=False, debug=False):
    """
    Perform comscripts unit testing.

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
        testdir = inspect.getfile(comscripts.tests)
        dirname = os.path.dirname(testdir)

        # Copy test data in "data" directory
        os.system('cp -r %s %s' % (dirname + '/data', 'data'))

        # Set 'TEST_DATA' and 'COMDATA' environment variable
        os.environ['TEST_DATA'] = 'data'
        os.environ['COMDATA']   = '%s/data/comptel/data' % (path)

    # ... otherwise set the calibration database to the one shipped with the
    # package; we don't need to set the 'TEST_DATA', this is done by the
    # test environment
    else:
        os.environ['CALDB']   = '%s/caldb'                  % (os.environ['TEST_SRCDIR'])
        os.environ['COMDATA'] = '%s/test/data/comptel/data' % (os.environ['TEST_SRCDIR'])

    # Create a local "pfiles" directory and set PFILES environment variable
    try:
        os.mkdir('pfiles')
    except:
        pass
    os.environ['PFILES'] = 'pfiles'

    # Copy the comscripts parameter files into the "pfiles" directory. For a
    # non-installed test we copy the parameter files from the "comscripts"
    # source directory into the "pfiles" directory, for an installed version
    # we get all parameter files from the "syspfiles" directory. Also make
    # sure that all parameter files are writable.
    if not installed:
        os.system('cp -r %s/modules/comscripts/*.par pfiles/' % (os.environ['TEST_SRCDIR']))
        os.system('chmod u+w pfiles/*')
    else:
        os.system('cp -r %s/syspfiles/*.par pfiles/' % (os.environ['CTOOLS']))
        os.system('chmod u+w pfiles/*')

    # Define list of test suites
    tests = [test_comobsselect.Test(),
             test_comobsbin.Test(),
             test_comobsback.Test(),
             test_comobsmodel.Test(),
             test_comobsadd.Test(),
             test_comobsres.Test(),
             test_comobssim.Test(),
             test_comlixfit.Test(),
             test_comlixmap.Test(),
             test_comsrcdetect.Test()]

    # Allocate test suite container
    suites = gammalib.GTestSuites('comscripts unit testing')

    # Set test suites and append them to suite container
    for suite in tests:
        suite.set()
        suites.append(suite)

    # Run test suite
    success = suites.run()

    # Save test results
    if not installed:
        suites.save('reports/comscripts.xml')
    else:
        suites.save('comscripts_reports.xml')

    # If debuging is requested then print test suites
    if debug:
        print(suites)

    # Remove temporary direction
    if installed:
        os.system('rm -rf %s' % (path))

    # Raise an exception in case of failure
    if not success:
        raise RuntimeError('At least one error occured during the test.')

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run tests
    test()
