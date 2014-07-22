#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctbkgcube tool.
#
# Copyright (C) 2014 Juergen Knoedlseder
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


# ============================= #
# Test class for ctbkgcube tool #
# ============================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctbkgcube tool.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set members
        self.events_name = "data/crab_events.fits"
        self.bkg_model   = "data/crab.xml"

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("ctbkgcube")

        # Append tests
        self.append(self.test_functional, "Test ctbkgcube functionality")

        # Return
        return

    # Test ctbkgcube functionnality
    def test_functional(self):
        """
        Test ctbkgcube functionnality.
        """
        # Set-up ctbkgcube
        bkgcube = ctools.ctbkgcube()
        bkgcube["infile"].filename(self.events_name)
        bkgcube["bkgmdl"].filename(self.bkg_model)
        bkgcube["cntmap"].filename("NONE")
        bkgcube["outfile"].filename("psfcube.fits")
        bkgcube["ebinalg"].string("LOG")
        bkgcube["emin"].real(0.1)
        bkgcube["emax"].real(100.0)
        bkgcube["enumbins"].integer(20)
        bkgcube["nxpix"].integer(10)
        bkgcube["nypix"].integer(10)
        bkgcube["binsz"].real(0.4)
        bkgcube["coordsys"].string("CEL")
        bkgcube["proj"].string("CAR")
        bkgcube["xref"].real(83.63)
        bkgcube["yref"].real(22.01)
        
        # Run tool
        self.test_try("Run ctbkgcube")
        try:
            bkgcube.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctbkgcube.")

        # Save background cube
        self.test_try("Save background cube")
        try:
            bkgcube.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving background cube.")

        # Return
        return
