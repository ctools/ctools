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

    # Allocate test suite container
    suites = gammalib.GTestSuites('ctools unit testing')

    # Allocate test suites and append them to the container
    suite_ctobssim    = test_ctobssim.Test()
    suite_ctselect    = test_ctselect.Test()
    suite_ctbin       = test_ctbin.Test()
    suite_ctlike      = test_ctlike.Test()
    suite_cttsmap     = test_cttsmap.Test()
    suite_ctmodel     = test_ctmodel.Test()
    suite_ctskymap    = test_ctskymap.Test()
    suite_ctexpcube   = test_ctexpcube.Test()
    suite_ctpsfcube   = test_ctpsfcube.Test()
    suite_ctedispcube = test_ctedispcube.Test()
    suite_ctbkgcube   = test_ctbkgcube.Test()
    suite_ctmapcube   = test_ctmapcube.Test()
    suite_ctcubemask  = test_ctcubemask.Test()
    suite_ctbutterfly = test_ctbutterfly.Test()
    suite_ctulimit    = test_ctulimit.Test()
    suite_cterror     = test_cterror.Test()
    suite_pipelines   = test_pipelines.Test()

    # Setup unit tests
    suite_ctobssim.set()
    suite_ctselect.set()
    suite_ctbin.set()
    suite_ctlike.set()
    suite_cttsmap.set()
    suite_ctmodel.set()
    suite_ctskymap.set()
    suite_ctexpcube.set()
    suite_ctpsfcube.set()
    suite_ctedispcube.set()
    suite_ctbkgcube.set()
    suite_ctmapcube.set()
    suite_ctcubemask.set()
    suite_ctbutterfly.set()
    suite_ctulimit.set()
    suite_cterror.set()
    suite_pipelines.set()

    # Append tests to container
    suites.append(suite_ctobssim)
    suites.append(suite_ctselect)
    suites.append(suite_ctbin)
    suites.append(suite_ctlike)
    suites.append(suite_cttsmap)
    suites.append(suite_ctmodel)
    suites.append(suite_ctskymap)
    suites.append(suite_ctexpcube)
    suites.append(suite_ctpsfcube)
    suites.append(suite_ctedispcube)
    suites.append(suite_ctbkgcube)
    suites.append(suite_ctmapcube)
    suites.append(suite_ctcubemask)
    suites.append(suite_ctbutterfly)
    suites.append(suite_ctulimit)
    suites.append(suite_cterror)
    suites.append(suite_pipelines)

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
