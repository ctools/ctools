#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctulimit tool.
#
# Copyright (C) 2015 Michael Mayer
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
# Test class for ctulimit tool #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctulimit tool.
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
        self.name("ctulimit")

        # Append tests
        self.append(self.test_functional, "Test ctulimit functionality")

        # Return
        return

    # Test cttsmap functionnality
    def test_functional(self):
        """
        Test ctulimit functionnality.
        """
        # Set-up cttsmap
        ulimit = ctools.ctulimit()
        ulimit["inobs"]   = self.events_name
        ulimit["inmodel"] = self.model_name
        ulimit["srcname"] = "Crab"
        ulimit["caldb"]   = self.caldb
        ulimit["irf"]     = self.irf
        
        # Run tool
        self.test_try("Run ctulimit")
        try:
            ulimit.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctulimit.")

        # Save results
        self.test_try("Save results")
        try:
            ulimit.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving results.")
        
        # Return
        return
