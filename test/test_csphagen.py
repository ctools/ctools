#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csphagen script.
#
# Copyright (C) 2017 Luigi Tibaldo
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
from testing import test


# ================================ #
# Test class for csphagen script #
# ================================ #
class Test(test):
    """
    Test class for csphagen script

    This test class makes unit tests for the csphagen script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csphagen')

        # Append tests
        self.append(self._test_cmd, 'Test csphagen on command line')
        self.append(self._test_python, 'Test csphagen from Python')

        # Return
        return

    # Test cslightcrv on command line
    def _test_cmd(self):
        """
        Test cslightcrv on the command line.
        """
        # Set script name
        csphagen = self._script('csphagen')

        cmd = csphagen +

