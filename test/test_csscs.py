#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csscs script.
#
# Copyright (C) 2020 Luigi Tibaldo
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
# Test class for csscs script #
# ========================= #
class Test(test):
    """
    Test class for csscs script

    This test class makes unit tests for the csscs script by using it from
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
        self.name('csscs')

        # Append tests
        self.append(self._test_cmd, 'Test csscs on command line')
        self.append(self._test_python, 'Test csscs from Python')
        self.append(self._test_pickeling, 'Test csscs pickeling')

        # Return
        return

    # Test csscs on command line
    def _test_cmd(self):
        """
        Test csscs on the command line
        """
        # Set script name
        csscs = self._script('csscs')

        # Setup csscs command
        cmd = csscs+' logfile="csscs_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Check csscs --help
        self._check_help(csscs)

        # Return
        return

    # Test csscs from Python
    def _test_python(self):
        """
        Test csscs from Python
        """
        # Allocate csscs
        script = cscripts.csscs()

        # Check that saving does not nothing
        script['logfile'] = 'csscs_py0.log'
        script.logFileOpen()
        script.save()
        self.test_assert(not os.path.isfile('csscs_py0.fits'),
             'Check that no FITS file has been created')

        # Check that clearing does not lead to an exception or segfault
        script.clear()

        # Now set csscs parameters
        script['logfile'] = 'csscs_py1.log'
        script['chatter'] = 2

        # Run csscs script
        script.logFileOpen()
        script.run()

        # Check result
        self._check_result(script)

        # Return
        return

    # Test csscs pickeling
    def _test_pickeling(self):
        """
        Test csscs pickeling
        """
        # Perform pickeling test of empty class
        self._pickeling(cscripts.csscs())

        # Setup script for pickling text
        script = cscripts.csscs()
        script['logfile'] = 'csscs_py1_pickle.log'
        script['chatter'] = 2

        # Perform pickeling tests of filled class
        obj = self._pickeling(script)

        # Run csscs script and save result
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # TODO: Check result file

        # Return
        return

    # Check csscs result
    def _check_result(self, script):
        """
        Check content of script

        Parameters
        ----------
        script : `~cscripts.csscs`
            csscs instance
        """
        # TODO: Implement test on script result

        # Return
        return
