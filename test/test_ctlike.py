#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctlike tool.
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


# ========================== #
# Test class for ctlike tool #
# ========================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctlike tool.
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
        self.name("ctlike")

        # Append tests
        self.append(self.test_functional, "Test ctlike functionality")

        # Return
        return

    # Test ctlike functionnality
    def test_functional(self):
        """
        Test ctlike functionnality.
        """
        # Set-up ctlike
        like = ctools.ctlike()
        like["inobs"]    = self.events_name
        like["inmodel"]  = self.model_name
        like["caldb"]    = self.caldb
        like["irf"]      = self.irf
        like["outmodel"] = "result.xml"

        # Run tool
        self.test_try("Run ctlike")
        try:
            like.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctlike.")

        # Save results
        self.test_try("Save results")
        try:
            like.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving results.")
        
        # Return
        return
