#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctskymap tool.
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


# ============================ #
# Test class for ctskymap tool #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctskymap tool.
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
        self.name("ctskymap")

        # Append tests
        self.append(self.test_functional, "Test ctskymap functionality")

        # Return
        return

    # Test ctskymap functionnality
    def test_functional(self):
        """
        Test ctskymap functionnality.
        """
        # Set-up ctskymap
        skymap = ctools.ctskymap()
        skymap["inobs"]    = self.events_name
        skymap["outmap"]   = "skymap.fits"
        skymap["emin"]     = 0.1
        skymap["emax"]     = 100
        skymap["nxpix"]    = 200
        skymap["nypix"]    = 200
        skymap["binsz"]    = 0.02
        skymap["coordsys"] = "CEL"
        skymap["proj"]     = "CAR"
        skymap["xref"]     = 83.63
        skymap["yref"]     = 22.01
        
        # Run tool
        self.test_try("Run ctskymap")
        try:
            skymap.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctskymap.")

        # Save counts cube
        self.test_try("Save sky map")
        try:
            skymap.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving sky map.")

        # Return
        return
