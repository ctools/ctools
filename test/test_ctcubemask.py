#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctcubemask tool.
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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


# ============================== #
# Test class for ctcubemask tool #
# ============================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctcubemask tool.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set members
        self.cntmap_name = "data/crab_cntmap.fits"
        self.regfile     = "data/exclusion.reg"

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("ctcubemask")

        # Append tests
        self.append(self.test_functional, "Test ctcubemask functionality")

        # Return
        return

    # Test ctcubemask functionnality
    def test_functional(self):
        """
        Test ctcubemask functionnality.
        """
        # Set-up ctcubemask
        mask = ctools.ctcubemask()
        mask["inobs"]   = self.cntmap_name
        mask["regfile"] = self.regfile
        mask["outcube"] = "filtered_cntmap.fits"
        mask["ra"]      = 83.63
        mask["dec"]     = 22.01
        mask["rad"]     = 2.0
        mask["emin"]    = 0.1
        mask["emax"]    = 100

        # Run tool
        self.test_try("Run ctcubemask")
        try:
            mask.run()
            self.test_try_success()
        except Exception, e:
            msg = "Exception occured in ctcubemask: %s." % (e,)
            self.test_try_failure(msg)

        # Save counts cube
        self.test_try("Save counts cube")
        try:
            mask.save()
            self.test_try_success()
        except Exception, e:
            msg = "Exception occured in saving counts cube: %s." % (e,)
            self.test_try_failure(msg)

        # Return
        return
