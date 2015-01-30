#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctselect tool.
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
# Test class for ctselect tool #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctselect tool.
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
        self.name("ctselect")

        # Append tests
        self.append(self.test_functional, "Test ctselect functionality")

        # Return
        return

    # Test ctselect functionnality
    def test_functional(self):
        """
        Test ctselect functionnality.
        """
        # Set-up ctselect
        select = ctools.ctselect()
        select["inobs"].filename(self.events_name)
        select["outobs"].filename("selected_events.fits")
        select["ra"].real(83.63)
        select["dec"].real(22.01)
        select["rad"].real(3.0)
        select["tmin"].real(0.0)
        select["tmax"].real(1800.0)
        select["emin"].real(0.1)
        select["emax"].real(100.0)

        # Run tool
        self.test_try("Run ctselect")
        try:
            select.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctselect.")

        # Save events
        self.test_try("Save events")
        try:
            select.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving events.")
        
        # Return
        return
