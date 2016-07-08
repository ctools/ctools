#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctlike tool.
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
from testing import test


# ========================== #
# Test class for ctlike tool #
# ========================== #
class Test(test):
    """
    Test class for ctlike tool
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
        self.name('ctlike')

        # Append tests
        self.append(self._test_cmd, 'Test ctlike on command line')
        self.append(self._test_python, 'Test ctlike from Python')

        # Return
        return

    # Test ctlike on command line
    def _test_cmd(self):
        """
        Test ctlike on the command line
        """
        # Set tool name
        ctlike = self._tool('ctlike')

        # Setup ctlike command
        cmd = ctlike+' inobs="data/crab_events.fits"'+ \
                     ' inmodel="data/crab.xml"'+ \
                     ' caldb="prod2" irf="South_0.5h"'+ \
                     ' outmodel="ctlike_cmd1.xml"'+ \
                     ' logfile="ctlike_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctlike_cmd1.xml')

        # Setup ctlike command
        cmd = ctlike+' inobs="event_file_that_does_not_exist.fits"'+ \
                     ' inmodel="data/crab.xml"'+ \
                     ' caldb="prod2" irf="South_0.5h"'+ \
                     ' outmodel="ctlike_cmd2.xml"'+ \
                     ' logfile="ctlike_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctlike from Python
    def _test_python(self):
        """
        Test ctlike from Python
        """
        # Set-up ctlike
        like = ctools.ctlike()
        like['inobs']    = 'data/crab_events.fits'
        like['inmodel']  = 'data/crab.xml'
        like['caldb']    = 'prod2'
        like['irf']      = 'South_0.5h'
        like['outmodel'] = 'ctlike_py1.xml'
        like['logfile']  = 'ctlike_py1.log'
        like['chatter']  = 2

        # Run ctlike tool
        like.logFileOpen()   # Make sure we get a log file
        like.run()
        like.save()

        # Check result file
        self._check_result_file('ctlike_py1.xml')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file
        """
        # Open result file
        result = gammalib.GModels(filename)

        # Check results
        self.test_value(result['Crab']['Prefactor'].value(), 1.58907e-16,
                        1.0e-3, 'Check fitted Crab Prefactor')
        self.test_value(result['Crab']['Prefactor'].error(), 0.0526982e-16,
                        1.0e-3, 'Check Crab Prefactor error')
        self.test_value(result['Crab']['Index'].value(), -2.43549,
                        1.0e-3, 'Check fitted Crab Index')
        self.test_value(result['Crab']['Index'].error(), 0.0248116,
                        1.0e-3, 'Check Crab Index error')
        self.test_value(result['Background']['Prefactor'].value(), 61.6919e-6,
                        1.0e-3, 'Check fitted background Prefactor')
        self.test_value(result['Background']['Prefactor'].error(), 1.49438e-6,
                        1.0e-3, 'Check background Prefactor error')
        self.test_value(result['Background']['Index'].value(), -2.20535,
                        1.0e-3, 'Check fitted background Index')
        self.test_value(result['Background']['Index'].error(), 0.0113269,
                        1.0e-3, 'Check background Index error')
        self.test_value(result['Background']['Sigma'].value(), 3.04252,
                        1.0e-3, 'Check fitted background Sigma')
        self.test_value(result['Background']['Sigma'].error(), 0.0307008,
                        1.0e-3, 'Check background Sigma error')

        # Return
        return
