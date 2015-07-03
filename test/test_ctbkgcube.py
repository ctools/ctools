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
        self.caldb      = "irf"
        self.irf        = "cta_dummy_irf"

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
        bkgcube["inobs"]    = self.events_name
        bkgcube["inmodel"]  = self.bkg_model
        bkgcube["incube"]   = "NONE"
        bkgcube["outcube"]  = "bkgcube.fits"
        bkgcube["outmodel"] = "bkgcube.xml"
        bkgcube["caldb"]    = self.caldb
        bkgcube["irf"]      = self.irf
        bkgcube["ebinalg"]  = "LOG"
        bkgcube["emin"]     = 0.1
        bkgcube["emax"]     = 100
        bkgcube["enumbins"] = 20
        bkgcube["nxpix"]    = 10
        bkgcube["nypix"]    = 10
        bkgcube["binsz"]    = 0.4
        bkgcube["coordsys"] = "CEL"
        bkgcube["proj"]     = "CAR"
        bkgcube["xref"]     = 83.63
        bkgcube["yref"]     = 22.01

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
