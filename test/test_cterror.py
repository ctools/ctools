#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cterror tool.
#
# Copyright (C) 2015-2016 Florent Forest
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


# =========================== #
# Test class for cterror tool #
# =========================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for cterror tool.
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
        self.model_name  = "data/crab.xml"
        self.result_name = "cterror_result.xml"
        self.caldb       = "irf"
        self.irf         = "cta_dummy_irf"

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('cterror')

        # Append tests
        self.append(self._test_cmd, 'Test cterror on command line')
        self.append(self._test_python, 'Test cterror from Python')

        # Return
        return

    # Test cterror on command line
    def _test_cmd(self):
        """
        Test cterror on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README.md'):
            cterror = '../src/cterror/cterror'
        else:
            cterror = 'cterror'

        # Setup cterror command
        cmd = cterror+' inobs="data/crab_events.fits"'+ \
                      ' inmodel="data/crab.xml" srcname="Crab"'+ \
                      ' caldb="prod2" irf="South_0.5h"'+ \
                      ' outmodel="cterror_cmd1.xml"'+ \
                      ' logfile="cterror_cmd1.log" chatter=1'

        # Execute cterror, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful cterror execution on command line')

        # Check result file
        self._check_result_file('cterror_cmd1.xml')

        # Setup cterror command
        cmd = cterror+' inobs="event_file_that_does_not_exist.fits"'+ \
                      ' inmodel="data/crab.xml" srcname="Crab"'+ \
                      ' caldb="prod2" irf="South_0.5h"'+ \
                      ' outmodel="cterror_cmd2.xml"'+ \
                      ' logfile="cterror_cmd2.log" chatter=1'

        # Execute cterror, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of cterror execution on command line')

        # Return
        return

    # Test cterror from Python
    def _test_python(self):
        """
        Test cterror from Python.
        """
        # Set-up cterror
        error = ctools.cterror()
        error['inobs']    = 'data/crab_events.fits'
        error['inmodel']  = 'data/crab.xml'
        error['outmodel'] = 'cterror_py1.xml'
        error['srcname']  = 'Crab'
        error['caldb']    = 'prod2'
        error['irf']      = 'South_0.5h'
        error['logfile']  = 'cterror_py1.log'
        error['chatter']  = 2

        # Run cterror tool
        error.logFileOpen()   # Make sure we get a log file
        error.run()
        error.save()

        # Check result file
        self._check_result_file('cterror_py1.xml')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file.
        """
        # Open result file
        result = gammalib.GModels(filename)

        # Check results
        self.test_value(result['Crab']['Prefactor'].value(), 1.58907e-16,
                        1.0e-3, 'Check fitted Crab Prefactor')
        self.test_value(result['Crab']['Prefactor'].error(), 0.0529105e-16,
                        1.0e-3, 'Check Crab Prefactor error')
        self.test_value(result['Crab']['Index'].value(), -2.43549,
                        1.0e-3, 'Check fitted Crab Index')
        self.test_value(result['Crab']['Index'].error(), 0.027804,
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
