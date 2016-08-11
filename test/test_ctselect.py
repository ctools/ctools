#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctselect tool.
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

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('ctselect')

        # Append tests
        self.append(self._test_cmd, 'Test ctselect on command line')
        self.append(self._test_python, 'Test ctselect from Python')

        # Return
        return

    # Test ctselect on command line
    def _test_cmd(self):
        """
        Test ctselect on the command line
        """
        # Set tool name
        ctselect = self._tool('ctselect')

        # Setup ctselect command
        cmd = ctselect+' inobs="'+self._events+'"'+ \
                       ' outobs="ctselect_cmd1.fits"'+ \
                       ' ra=83.63 dec=22.01 rad=3.0'+ \
                       ' tmin=0.0 tmax=1800.0'+ \
                       ' emin=0.1 emax=100.0'+ \
                       ' logfile="ctselect_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctselect_cmd1.fits')

        # Setup ctselect command
        cmd = ctselect+' inobs="event_file_that_does_not_exist.fits"'+ \
                       ' outobs="ctselect_cmd2.fits"'+ \
                       ' ra=83.63 dec=22.01 rad=3.0'+ \
                       ' tmin=0.0 tmax=1800.0'+ \
                       ' emin=0.1 emax=100.0'+ \
                       ' logfile="ctselect_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctselect from Python
    def _test_python(self):
        """
        Test ctselect from Python
        """
        # Allocate ctselect
        select = ctools.ctselect()

        # Check that empty ctselect tool holds an empty observation
        self._check_obs(select.obs(), nobs=0)

        # Check that saving does nothing
        select['outobs']  = 'ctselect_py0.fits'
        select['logfile'] = 'ctselect_py0.log'
        select.logFileOpen()
        select.save()
        self.test_assert(not os.path.isfile('ctselect_py0.fits'),
             'Check that no event list has been created')

        # Check that clearing does not lead to an exception or segfault
        select.clear()

        # Now set ctselect parameters
        select = ctools.ctselect()
        select['inobs']   = self._events
        select['ra']      = 83.63
        select['dec']     = 22.01
        select['rad']     = 3
        select['tmin']    = 0
        select['tmax']    = 1800
        select['emin']    = 0.1
        select['emax']    = 100
        select['outobs']  = 'ctselect_py1.fits'
        select['logfile'] = 'ctselect_py1.log'
        select['chatter'] = 2

        # Run ctselect tool
        select.logFileOpen()   # Make sure we get a log file
        select.run()
        select.save()

        # Check result file
        self._check_result_file('ctselect_py1.fits')

        # Copy ctselect tool
        cpy_select = select.copy()

        # Check observation of ctselect copy
        self._check_obs(cpy_select.obs())

        # Execute copy of ctselect tool again, now with a higher chatter
        # level than before and without any selection
        cpy_select['ra']      = 'NONE'
        cpy_select['dec']     = 'NONE'
        cpy_select['rad']     = 'NONE'
        cpy_select['tmin']    = 'NONE'
        cpy_select['tmax']    = 'NONE'
        cpy_select['emin']    = 'NONE'
        cpy_select['emax']    = 'NONE'
        cpy_select['outobs']  = 'ctselect_py2.fits'
        cpy_select['logfile'] = 'ctselect_py2.log'
        cpy_select['chatter'] = 3
        cpy_select.logFileOpen()  # Needed to get a new log file
        cpy_select.execute()

        # Check result file
        self._check_result_file('ctselect_py2.fits')
        print('Step3')

        # Execute again the copy of ctselect tool again, now using the pointing
        # information in the input observation
        cpy_select['usepnt']  = True
        cpy_select['outobs']  = 'ctselect_py3.fits'
        cpy_select['logfile'] = 'ctselect_py3.log'
        cpy_select['chatter'] = 4
        cpy_select.logFileOpen()  # Needed to get a new log file
        cpy_select.execute()

        # Check result file
        self._check_result_file('ctselect_py3.fits')

        # Now clear copy of ctselect tool
        cpy_select.clear()

        # Check that cleared ctselect tool holds no observations and events
        self._check_obs(cpy_select.obs(), nobs=0)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file
        """
        # Open result file
        events = gammalib.GCTAEventList(filename)

        # Check event list
        self._check_events(events, nevents=6127)

        # Return
        return

    # Check observation and event list
    def _check_obs(self, obs, nobs=1, nevents=6127):
        """
        Check observation and event list

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Models
        nobs : int, optional
            Expected number of observations
        nevents : int, optional
            Expected number of events
        """
        # Check number of observations
        self.test_value(obs.size(), nobs, 'Check number of observations')

        # If there is an observation then check the event list of the first
        # one
        if obs.size() > 0:
            self._check_events(obs[0].events(), nevents=nevents)

        # Return
        return

    # Check events
    def _check_events(self, events, nevents=6127):
        """
        Check event list

        Parameters
        ----------
        events : `~gammalib.GCTAEventList`
            Event list
        nevents : int, optional
            Expected number of events
        """
        # Check event list
        self.test_value(events.size(), nevents, 'Check number of events')
        self.test_value(events.roi().centre().dir().ra_deg(), 83.63, 1.0e-6,
                        'Check for RoI centre Right Ascension')
        self.test_value(events.roi().centre().dir().dec_deg(), 22.01, 1.0e-6,
                        'Check for RoI centre Declination')
        self.test_value(events.roi().radius(), 3.0, 1.0e-6,
                        'Check for RoI radius')

        # Return
        return
