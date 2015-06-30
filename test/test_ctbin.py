#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctbin tool.
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


# ========================= #
# Test class for ctbin tool #
# ========================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctbin tool.
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

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("ctbin")

        # Append tests
        self.append(self.test_functional, "Test ctbin functionality")

        # Return
        return

    # Test ctbin functionnality
    def test_functional(self):
        """
        Test ctbin functionnality.
        """
        # Set-up ctbin
        bin = ctools.ctbin()
        bin["inobs"]    = self.events_name
        bin["outcube"]  = "cntmap.fits"
        bin["ebinalg"]  = "LOG"
        bin["emin"]     = 0.1
        bin["emax"]     = 100.0
        bin["enumbins"] = 20
        bin["nxpix"]    = 200
        bin["nypix"]    = 200
        bin["binsz"]    = 0.02
        bin["coordsys"] = "CEL"
        bin["proj"]     = "CAR"
        bin["xref"]     = 83.63
        bin["yref"]     = 22.01
        
        # Run tool
        self.test_try("Run ctbin")
        try:
            bin.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctbin.")

        # Save counts cube
        self.test_try("Save counts cube")
        try:
            bin.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving counts cube.")

        # Return
        return
