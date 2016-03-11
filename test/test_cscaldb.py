#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cscaldb script.
#
# Copyright (C) 2016 Juergen Knoedlseder
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
import cscripts


# ============================= #
# Test class for cscaldb script #
# ============================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for cscaldb script.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("cscaldb")

        # Append tests
        self.append(self._test_functional, "Test cscaldb functionality")

        # Return
        return

    # Test cscaldb functionality
    def _test_functional(self):
        """
        Test cscaldb functionality.
        """
        # Set-up cscaldb
        caldb = cscripts.cscaldb()

        # Run script
        self.test_try("Run cscaldb")
        try:
            caldb.run()
            self.test_try_success()
        except:
            msg = "Exception occured in cscaldb."
            self.test_try_failure(msg)

        # Return
        return
