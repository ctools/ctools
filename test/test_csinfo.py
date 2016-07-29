#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csinfo script.
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
from testing import test


# ============================ #
# Test class for csinfo script #
# ============================ #
class Test(test):
    """
    Test class for csinfo script
    """

    # Constructor
    def __init__(self):
        """
        Constructor.
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
        self.name('csinfo')

        # Append tests
        self.append(self._test_cmd, 'Test csinfo on command line')

        # Return
        return

    # Test csinfo on command line
    def _test_cmd(self):
        """
        Test csinfo on the command line.
        """
        # Set script name
        csinfo = self._script('csinfo')

        # Setup csinfo command
        cmd = csinfo

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Execute script without arguments
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Execute script with "list" argument
        self.test_assert(self._execute(cmd+' list') == 0,
             'Check successful execution with "list" argument from command line')

        # Execute script with "check" argument
        self.test_assert(self._execute(cmd+' check') == 0,
             'Check successful execution with "check" argument from command line')

        # Execute script with "info" argument
        self.test_assert(self._execute(cmd+' info') == 0,
             'Check successful execution with "info" argument from command line')

        # Execute script with invalid argument
        self.test_assert(self._execute(cmd+' not_a_valid_argument') != 0,
             'Check unsuccessful execution with unsupported argument from '
             'command line')

        # Return
        return
