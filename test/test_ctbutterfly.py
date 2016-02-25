#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctbutterfly tool.
#
# Copyright (C) 2014-2016 Michal Mayer
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


# =============================== #
# Test class for ctbutterfly tool #
# =============================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctbutterfly tool.
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
        self.name("ctbutterfly")

        # Append tests
        self.append(self.test_functional, "Test ctbutterfly functionality")

        # Return
        return

    # Test ctbutterfly functionnality
    def test_functional(self):
        """
        Test ctbutterfly functionnality.
        """
        # Set-up ctbutterfly
        butterfly = ctools.ctbutterfly()
        butterfly["inobs"]   = self.events_name
        butterfly["inmodel"] = self.model_name
        butterfly["srcname"] = "Crab"
        butterfly["caldb"]   = self.caldb
        butterfly["irf"]     = self.irf
        butterfly["emin"]    = 0.1
        butterfly["emax"]    = 100
        butterfly["outfile"] = "butterfly.txt"

        # Run tool
        self.test_try("Run ctbutterfly")
        try:
            butterfly.run()
            self.test_try_success()
        except Exception as e:
            msg = "Exception occured in ctbutterfly: %s." % (e,)
            self.test_try_failure(msg)

        # Save results
        self.test_try("Save results")
        try:
            butterfly.save()
            self.test_try_success()
        except Exception as e:
            msg = "Exception occured in saving butterfly diagram: %s." % (e,)
            self.test_try_failure(msg)

        # Return
        return
