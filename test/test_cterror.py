#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cterror tool.
#
# Copyright (C) 2015 Florent Forest
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


# =========================== #
# Test class for cterror tool #
# =========================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for cterror tool.
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
        self.model_name  = "data/crab.xml"
        self.result_name = "cterror_result.xml"
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
        self.name("cterror")

        # Append tests
        self.append(self.test_functional, "Test cterror functionality")

        # Return
        return

    # Test cttsmap functionnality
    def test_functional(self):
        """
        Test cterror functionnality.
        """
        # Set-up cterror
        error = ctools.cterror()
        error["inobs"].filename(self.events_name)
        error["inmodel"].filename(self.model_name)
        error["outmodel"].filename(self.result_name)
        error["srcname"].string("Crab")
        error["caldb"].string(self.caldb)
        error["irf"].string(self.irf)
        
        # Run tool
        self.test_try("Run cterror")
        try:
            error.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in cterror.")

        # Save results
        self.test_try("Save results")
        try:
            error.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving results.")
        
        # Return
        return
