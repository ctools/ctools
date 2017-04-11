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
        self._phase_file     = self._datadir + '/model_temporal_phasecurve.fits'
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
                      ' inmodel="'+self._model_file+'" source="Crab"'+ \
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
                      ' inmodel="'+self._model_file+'" source="Crab"'+ \
                      ' logfile="ctphase_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctphase --help option
        cmd = ctphase+' --help'

        # Check if execution was successful in case that the CTOOLS
        # environment variable was set or failed otherwise
        if 'CTOOLS' in os.environ:
            self.test_value(self._execute(cmd), 0,
                 'Check successful execution with --help option')
        else:
            self.test_assert(self._execute(cmd) != 0,
                 'Check execution failure with --help option')

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
        phase['source']  = 'Crab'
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
        phase['source']  = 'Crab'
        phase['outobs']  = 'ctphase_py1.fits'
        phase['logfile'] = 'ctphase_py1.log'
        phase['chatter'] = 2

        # Run ctphase tool
        phase.logFileOpen()   # Make sure we get a log file
        phase.run()
        phase.save()

        # Check result file
        self._check_result_file('ctphase_py1.fits')

        # Now run ctphase tool without a model file
        phase['inmodel']    = 'NONE'
        phase['source']     = 'NONE'
        phase['phasecurve'] = self._phase_file
        phase['mjd']        = 51544.5
        phase['phase']      = 0.0
        phase['f0']         = 1.0
        phase['f1']         = 0.1
        phase['f2']         = 0.01
        phase['outobs']     = 'ctphase_py2.fits'
        phase['logfile']    = 'ctphase_py2.log'
        phase['chatter']    = 3

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

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nevents=6141):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Event list file name
        nevents : int, optional
            Expected number of events
        """
        # Open result file
        events      = gammalib.GCTAEventList(filename)
        phases      = [events[i].phase() for i in range(events.size())]
        test_result = int(sum(phases)>0)
        self.test_value(test_result, 1, 'Check number of events')

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
