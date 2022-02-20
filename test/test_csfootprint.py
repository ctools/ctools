#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csfootprint script.
#
# Copyright (C) 2022 Juergen Knoedlseder
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


# ================================= #
# Test class for csfootprint script #
# ================================= #
class Test(test):
    """
    Test class for csfootprint script

    This test class makes unit tests for the csfootprint script by using it from
    the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set members
        self._infile = self._datadir + '/statistics.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csfootprint')

        # Append tests
        self.append(self._test_cmd, 'Test csfootprint on command line')
        self.append(self._test_python, 'Test csfootprint from Python')
        self.append(self._test_pickeling, 'Test csfootprint pickeling')

        # Return
        return

    # Test csfootprint on command line
    def _test_cmd(self):
        """
        Test csfootprint on the command line
        """
        # Set script name
        csfootprint = self._script('csfootprint')

        # Setup csfootprint command
        cmd = csfootprint+' infile="'+self._infile+'" outfile=NONE tmin=NONE'+ \
                          ' tmax=NONE logfile="csfootprint_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Setup csfootprint command
        cmd = csfootprint+' infile="file_that_does_not_exist.xml" outfile=NONE tmin=NONE'+ \
                          ' tmax=NONE logfile="csfootprint_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check csfootprint --help
        self._check_help(csfootprint)

        # Return
        return

    # Test csfootprint from Python
    def _test_python(self):
        """
        Test csfootprint from Python
        """
        # Allocate csfootprint
        script = cscripts.csfootprint()

        # Check that saving does not nothing
        script['infile']  = self._infile
        script['outfile'] = 'NONE'
        script['tmin']    = 'NONE'
        script['tmax']    = 'NONE'
        script['logfile'] = 'csfootprint_py0.log'
        script['chatter'] = 2

        # Run csfootprint script
        script.logFileOpen()
        script.run()

        # Check result
        self._check_result(script)

        # Return
        return

    # Test csfootprint pickeling
    def _test_pickeling(self):
        """
        Test csfootprint pickeling
        """
        # Perform pickeling test of empty class
        self._pickeling(cscripts.csfootprint())

        # Setup script for pickling text
        script = cscripts.csfootprint()
        script['infile']  = self._infile
        script['outfile'] = 'NONE'
        script['tmin']    = 'NONE'
        script['tmax']    = 'NONE'
        script['logfile'] = 'csfootprint_py1_pickle.log'
        script['chatter'] = 2

        # Perform pickeling tests of filled class
        obj = self._pickeling(script)

        # Run csfootprint script and save result
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # Check result
        self._check_result(script)

        # Return
        return

    # Check csfootprint result
    def _check_result(self, script):
        """
        Check content of script

        Parameters
        ----------
        script : `~cscripts.csfootprint`
            csfootprint instance
        """
        # TODO: Implement test on script result

        # Return
        return
