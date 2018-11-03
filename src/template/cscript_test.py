#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the xxx script.
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
import cscripts
from testing import test


# ========================= #
# Test class for xxx script #
# ========================= #
class Test(test):
    """
    Test class for xxx script

    This test class makes unit tests for the xxx script by using it from
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
        self.append(self._test_pickeling, 'Test xxx pickeling')

        # Return
        return

    # Test xxx on command line
    def _test_cmd(self):
        """
        Test xxx on the command line
        """
        # Set script name
        xxx = self._script('xxx')

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
        script = cscripts.xxx()

        # Check that saving does not nothing
        script['logfile'] = 'xxx_py0.log'
        script.logFileOpen()
        script.save()
        self.test_assert(not os.path.isfile('xxx_py0.fits'),
             'Check that no FITS file has been created')

        # Check that clearing does not lead to an exception or segfault
        script.clear()

        # Now set xxx parameters
        script['logfile'] = 'xxx_py1.log'
        script['chatter'] = 2

        # Run xxx script
        script.logFileOpen()
        script.run()

        # Check result
        self._check_result(script)

        # Return
        return

    # Test xxx pickeling
    def _test_pickeling(self):
        """
        Test xxx pickeling
        """
        # Perform pickeling test of empty class
        self._pickeling(cscripts.xxx())

        # Setup script for pickling text
        script = cscripts.xxx()
        script['logfile'] = 'xxx_py1_pickle.log'
        script['chatter'] = 2

        # Perform pickeling tests of filled class
        obj = self._pickeling(script)

        # Run xxx script and save result
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # TODO: Check result file

        # Return
        return

    # Check xxx result
    def _check_result(self, script):
        """
        Check content of script

        Parameters
        ----------
        script : `~cscripts.xxx`
            xxx instance
        """
        # TODO: Implement test on script result

        # Return
        return
