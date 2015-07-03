#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctpsfcube tool.
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
# Test class for ctpsfcube tool #
# ============================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctpsfcube tool.
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
        self.caldb       = "irf"
        self.irf         = "cta_dummy_irf"

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("ctpsfcube")

        # Append tests
        self.append(self.test_functional, "Test ctpsfcube functionality")

        # Return
        return

    # Test ctpsfcube functionnality
    def test_functional(self):
        """
        Test ctpsfcube functionnality.
        """
        # Set-up ctpsfcube
        psfcube = ctools.ctpsfcube()
        psfcube["inobs"]    = self.events_name
        psfcube["incube"]   = "NONE"
        psfcube["outcube"]  = "psfcube.fits"
        psfcube["caldb"]    = self.caldb
        psfcube["irf"]      = self.irf
        psfcube["ebinalg"]  = "LOG"
        psfcube["emin"]     = 0.1
        psfcube["emax"]     = 100
        psfcube["enumbins"] = 20
        psfcube["nxpix"]    = 10
        psfcube["nypix"]    = 10
        psfcube["binsz"]    = 0.4
        psfcube["coordsys"] = "CEL"
        psfcube["proj"]     = "CAR"
        psfcube["xref"]     = 83.63
        psfcube["yref"]     = 22.01
        psfcube["amax"]     = 0.3
        psfcube["anumbins"] = 10

        # Run tool
        self.test_try("Run ctpsfcube")
        try:
            psfcube.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctpsfcube.")

        # Save PSF cube
        self.test_try("Save PSF cube")
        try:
            psfcube.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving PSF cube.")

        # Return
        return
