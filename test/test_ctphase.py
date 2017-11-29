#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctphase tool.
#
# Copyright (C) 2017 by Joshua Cardenzana
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


# ============================ #
# Test class for ctselect tool #
# ============================ #
class Test(test):
    """
    Test class for ctselect tool
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test data
        self._model_file     = self._datadir + '/model_temporal_phasecurve.xml'
        self._invalid_events = self._datadir + '/invalid_event_list.fits'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('ctphase')

        # Append tests
        self.append(self._test_cmd, 'Test ctphase on command line')
        self.append(self._test_python, 'Test ctphase from Python')

        # Return
        return

    # Test ctphase on command line
    def _test_cmd(self):
        """
        Test ctphase on the command line
        """
        # Set tool name
        ctphase = self._tool('ctphase')

        # Setup ctphase command
        cmd = ctphase+' inobs="'+self._events+'"'+ \
                      ' outobs="ctphase_cmd1.fits"'+ \
                      ' inmodel="'+self._model_file+'" srcname="Crab"'+ \
                      ' logfile="ctphase_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctphase_cmd1.fits')

        # Setup ctphase command
        cmd = ctphase+' inobs="event_file_that_does_not_exist.fits"'+ \
                      ' outobs="ctphase_cmd2.fits"'+ \
                      ' inmodel="'+self._model_file+'" srcname="Crab"'+ \
                      ' logfile="ctphase_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Check ctphase --help
        self._check_help(ctphase)

        # Return
        return

    # Test ctphase from Python
    def _test_python(self):
        """
        Test ctphase from Python
        """
        # Allocate empty ctphase tool
        phase = ctools.ctphase()

        # Check that empty ctphase tool holds an empty observation
        self._check_obs(phase.obs(), nobs=0)

        # Check that saving does nothing
        phase['inmodel'] = self._model_file
        phase['srcname'] = 'Crab'
        phase['outobs']  = 'ctphase_py0.fits'
        phase['logfile'] = 'ctphase_py0.log'
        phase.logFileOpen()
        phase.save()
        self.test_assert(not os.path.isfile('ctphase_py0.fits'),
             'Check that no event list has been created')

        # Check that clearing does not lead to an exception or segfault
        phase.clear()

        # Now set ctphase parameters
        phase['inobs']   = self._events
        phase['inmodel'] = self._model_file
        phase['srcname'] = 'Crab'
        phase['outobs']  = 'ctphase_py1.fits'
        phase['logfile'] = 'ctphase_py1.log'
        phase['chatter'] = 2

        # Run ctphase tool
        phase.logFileOpen()   # Make sure we get a log file
        phase.run()
        phase.save()

        # Check result file
        self._check_result_file('ctphase_py1.fits')

        # Clear phase tool
        phase.clear()

        # Now run ctphase tool without a model file
        phase['inobs']   = self._events
        phase['inmodel'] = 'NONE'
        phase['srcname'] = 'NONE'
        phase['mjd']     = 58849.0
        phase['phase']   = 0.0
        phase['f0']      = 1.0
        phase['f1']      = 0.1
        phase['f2']      = 0.01
        phase['outobs']  = 'ctphase_py2.fits'
        phase['logfile'] = 'ctphase_py2.log'
        phase['chatter'] = 3

        # Execute ctphase tool
        phase.logFileOpen()   # Make sure we get a log file
        phase.execute()

        # Check result file
        self._check_result_file('ctphase_py2.fits')

        # Copy ctphase tool
        cpy_phase = phase.copy()

        # Clear phase tool
        phase.clear()

        # Execute copy
        cpy_phase['outobs']  = 'ctphase_py3.fits'
        cpy_phase['logfile'] = 'ctphase_py3.log'
        cpy_phase['chatter'] = 4
        cpy_phase['publish'] = True
        cpy_phase.logFileOpen()  # Needed to get a new log file
        cpy_phase.execute()

        # Check result file
        self._check_result_file('ctphase_py3.fits')

        # Now run ctphase tool without a model file and a significant
        # MJD offset
        phase['inobs']   = self._events
        phase['inmodel'] = 'NONE'
        phase['srcname'] = 'NONE'
        phase['mjd']     = 58860.0
        phase['phase']   = 0.0
        phase['f0']      = 1.0
        phase['f1']      = 0.1
        phase['f2']      = 0.0001
        phase['outobs']  = 'ctphase_py4.fits'
        phase['logfile'] = 'ctphase_py4.log'
        phase['chatter'] = 3

        # Execute ctphase tool
        phase.logFileOpen()   # Make sure we get a log file
        phase.execute()

        # Check result file
        self._check_result_file('ctphase_py4.fits')

        # Test invalid source name
        self.test_try('Test ctphase with invalid source name')
        try:
            phase = ctools.ctphase()
            phase['inobs']   = self._events
            phase['inmodel'] = self._model_file
            phase['srcname'] = 'Venus'
            phase['outobs']  = 'ctphase_py5.fits'
            phase['logfile'] = 'ctphase_py5.log'
            phase.logFileOpen()
            phase.execute()
            self.test_try_failure('Exception not thrown')
        except (ValueError):
            self.test_try_success()

        # Test invalid model type
        self.test_try('Test ctphase with invalid model')
        try:
            phase = ctools.ctphase()
            phase['inobs']   = self._events
            phase['inmodel'] = self._model
            phase['srcname'] = 'Background'
            phase['outobs']  = 'ctphase_py6.fits'
            phase['logfile'] = 'ctphase_py6.log'
            phase.logFileOpen()
            phase.execute()
            self.test_try_failure('Exception not thrown')
        except (ValueError):
            self.test_try_success()

        # Test invalid model without temporal phase curve model
        self.test_try('Test ctphase with invalid model')
        try:
            phase = ctools.ctphase()
            phase['inobs']   = self._events
            phase['inmodel'] = self._model
            phase['srcname'] = 'Crab'
            phase['outobs']  = 'ctphase_py7.fits'
            phase['logfile'] = 'ctphase_py7.log'
            phase.logFileOpen()
            phase.execute()
            self.test_try_failure('Exception not thrown')
        except (ValueError):
            self.test_try_success()

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nevents=22220):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Event list file name
        nevents : int, optional
            Expected number of events
        """
        # Open event list
        events = gammalib.GCTAEventList(filename)

        # Check event list
        self.test_value(events.size(), nevents, 'Check number of events')

        # Check phases
        phases      = [events[i].phase() for i in range(events.size())]
        test_result = int(sum(phases) > 0)
        self.test_value(test_result, 1, 'Check whether phase information is set')

        # Return
        return

    # Check observation and event list
    def _check_obs(self, obs, nobs=1):
        """
        Check observation and event list

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Models
        nobs : int, optional
            Expected number of observations
        """
        # Check number of observations
        self.test_value(obs.size(), nobs, 'Check number of observations')

        # Return
        return
