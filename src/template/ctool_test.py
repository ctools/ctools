#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the xxx tool.
#
# Copyright (C) [YEAR] [AUTHOR]
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
import os
import gammalib
import ctools
from testing import test


# ======================= #
# Test class for xxx tool #
# ======================= #
class Test(test):
    """
    Test class for xxx tool

    This test class makes unit tests for the xxx tool by using it from
    the command line and from Python.
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
        Set all test functions
        """
        # Set test name
        self.name('xxx')

        # Append tests
        self.append(self._test_cmd, 'Test xxx on command line')
        self.append(self._test_python, 'Test xxx from Python')

        # Return
        return

    # Test xxx on command line
    def _test_cmd(self):
        """
        Test xxx on the command line
        """
        # Set tool name
        xxx = self._tool('xxx')

        # Setup xxx command
        cmd = xxx+' logfile="xxx_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Check xxx --help
        self._check_help(xxx)

        # Return
        return

    # Test xxx from Python
    def _test_python(self):
        """
        Test xxx from Python
        """
        # Allocate xxx
        tool = ctools.xxx()

        # Check that saving does not nothing
        tool['logfile'] = 'xxx_py0.log'
        tool.logFileOpen()
        tool.save()
        self.test_assert(not os.path.isfile('xxx_py0.fits'),
             'Check that no FITS file has been created')

        # Check that clearing does not lead to an exception or segfault
        tool.clear()

        # Now set xxx parameters
        tool['logfile'] = 'xxx_py1.log'
        tool['chatter'] = 2

        # Run xxx tool
        tool.logFileOpen()
        tool.run()

        # Check result
        self._check_result(tool)

        # Copy xxx tool
        cpy_tool = tool.copy()

        # Check result
        self._check_result(cpy_tool)

        # Run copy of xxx tool again
        cpy_tool['logfile'] = 'xxx_py2.log'
        cpy_tool['chatter'] = 3
        cpy_tool.logFileOpen()
        cpy_tool.run()

        # TODO: Check result

        # Return
        return

    # Check xxx result
    def _check_result(self, tool):
        """
        Check content of an observation

        Parameters
        ----------
        tool : `~ctools.xxx`
            xxx instance
        """
        # TODO: Implement test on tool result

        # Return
        return
