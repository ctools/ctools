#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctfindvar tool.
#
# Copyright (C) 2018 Simon Bonnefoy
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
# Test class for ctfindvar tool #
# ======================= #
class Test(test):
    """
    Test class for ctfindvar tool

    This test class makes unit tests for the ctfindvar tool by using it from
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
        self.name('ctfindvar')

        # Append tests
        self.append(self._test_cmd, 'Test ctfindvar on command line')
        self.append(self._test_python, 'Test ctfindvar from Python')

        # Return
        return

    # Test ctfindvar on command line
    def _test_cmd(self):
        """
        Test ctfindvar on the command line
        """
        # Set tool name
        ctfindvar = self._tool('ctfindvar')

        # Setup ctfindvar command
        cmd = ctfindvar+' logfile="ctfindvar_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Check ctfindvar --help
        self._check_help(ctfindvar)

        # Return
        return

    # Test ctfindvar from Python
    def _test_python(self):
        """
        Test ctfindvar from Python
        """
        # Allocate ctfindvar
        tool = ctools.ctfindvar()

        # Check that saving does not nothing
        tool['logfile'] = 'ctfindvar_py0.log'
        tool.logFileOpen()
        tool.save()
        self.test_assert(not os.path.isfile('ctfindvar_py0.fits'),
             'Check that no FITS file has been created')

        # Check that clearing does not lead to an exception or segfault
        tool.clear()

        # Now set ctfindvar parameters
        tool['logfile'] = 'ctfindvar_py1.log'
        tool['chatter'] = 2

        # Run ctfindvar tool
        tool.logFileOpen()
        tool.run()

        # Check result
        self._check_result(tool)

        # Copy ctfindvar tool
        cpy_tool = tool.copy()

        # Check result
        self._check_result(cpy_tool)

        # Run copy of ctfindvar tool again
        cpy_tool['logfile'] = 'ctfindvar_py2.log'
        cpy_tool['chatter'] = 3
        cpy_tool.logFileOpen()
        cpy_tool.run()

        # Check result
        self._check_result(cpy_tool)

        # Return
        return

    # Check ctfindvar result
    def _check_result(self, tool):
        """
        Check content of tool

        Parameters
        ----------
        tool : `~ctools.ctfindvar`
            ctfindvar instance
        """
        # TODO: Implement test on tool result

        # Return
        return
