#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctools package.
#
# Copyright (C) 2012-2014 Juergen Knoedlseder
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
import test_ctobssim
import test_ctselect
import test_ctbin
import test_ctlike
import test_cttsmap
import test_ctmodel
import test_ctskymap
import test_ctexpcube
import test_ctpsfcube
import test_ctbkgcube
import test_ctcubemask
import test_ctbutterfly
import test_pipelines


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Perform unit testing for Python interface.
    """
    # Allocate test suite container
    suites = gammalib.GTestSuites("ctools unit testing")

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
    suite_ctbkgcube   = test_ctbkgcube.Test()
    suite_ctcubemask  = test_ctcubemask.Test()
    suite_ctbutterfly = test_ctbutterfly.Test()
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
    suite_ctbkgcube.set()
    suite_ctcubemask.set()
    suite_ctbutterfly.set()
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
    suites.append(suite_ctbkgcube)
    suites.append(suite_ctcubemask)
    suites.append(suite_ctbutterfly)
    suites.append(suite_pipelines)

    # Set PFILES environment variable
    try:
        os.mkdir("pfiles")
    except:
        pass
    os.system("cp -r ../src/*/*.par pfiles/")
    os.environ['PFILES'] = "pfiles"

    # Run test suite
    success = suites.run()

    # Save test results
    suites.save("reports/ctools.xml")

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Exit with return code
    sys.exit(rc)
