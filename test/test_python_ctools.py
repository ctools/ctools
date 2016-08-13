#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctools package
#
# Copyright (C) 2012-2016 Juergen Knoedlseder
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
import ctools
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import test_ctobssim
import test_ctselect
import test_ctbin
import test_ctlike
import test_cttsmap
import test_ctmodel
import test_ctskymap
import test_ctexpcube
import test_ctpsfcube
import test_ctedispcube
import test_ctbkgcube
import test_ctmapcube
import test_ctcubemask
import test_ctbutterfly
import test_ctulimit
import test_cterror
import test_pipelines


# ========================= #
# Perform ctools unit tests #
# ========================= #
def test(installed=False):
    """
    Perform unit testing for ctools.

    Parameters
    ----------
    installed : bool, optional
        Flag indicating whether the script has been installed or not
    """
    # If we have an installed version then create a temporary
    # directory and copy over all information that is needed
    if installed:

        # Create temporary working directory
        import tempfile
        path = tempfile.mkdtemp()
        os.chdir(path)

        # Get test directory
        import inspect
        testdir = inspect.getfile(ctools.tests)
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

    # Copy the ctools parameter files into the "pfiles" directory. For a
    # non-installed test we copy the parameter files from the respective
    # source directories into the "pfiles" directory, for an installed version
    # we get all parameter files from the "syspfiles" directory. Also make
    # sure that all parameter files are writable.
    if not installed:
        os.system('cp -r %s/src/*/*.par pfiles/' % (os.environ['TEST_SRCDIR']))
        os.system('chmod u+w pfiles/*')
    else:
        os.system('cp -r %s/syspfiles/*.par pfiles/' % (os.environ['CTOOLS']))
        os.system('chmod u+w pfiles/*')

    # Define list of test suites
    tests = [test_ctobssim.Test(),
             test_ctselect.Test(),
             test_ctbin.Test(),
             test_ctlike.Test(),
             test_cttsmap.Test(),
             test_ctmodel.Test(),
             test_ctskymap.Test(),
             test_ctexpcube.Test(),
             test_ctpsfcube.Test(),
             test_ctedispcube.Test(),
             test_ctbkgcube.Test(),
             test_ctmapcube.Test(),
             test_ctcubemask.Test(),
             test_ctbutterfly.Test(),
             test_ctulimit.Test(),
             test_cterror.Test(),
             test_pipelines.Test()]

    # Allocate test suite container
    suites = gammalib.GTestSuites('ctools unit testing')

    # Set test suites and append them to suite container
    for suite in tests:
        suite.set()
        suites.append(suite)

    # Run test suite
    success = suites.run()

    # If we have a non-installed version then save test results
    if not installed:
        suites.save('reports/ctools.xml')

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
