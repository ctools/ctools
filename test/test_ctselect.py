#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctselect tool.
#
# Copyright (C) 2014-2018 Juergen Knoedlseder
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
        self._phased_events  = self._datadir + '/phased_events.fits'
        self._invalid_events = self._datadir + '/invalid_event_list.fits'

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
                       ' rad=2.0'+ \
                       ' tmin=2020-01-01T00:01:10 tmax=2020-01-01T00:03:40'+ \
                       ' emin=0.2 emax=80.0'+ \
                       ' logfile="ctselect_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctselect_cmd1.fits')

        # Setup ctselect command
        cmd = ctselect+' inobs="event_file_that_does_not_exist.fits"'+ \
                       ' outobs="ctselect_cmd2.fits"'+ \
                       ' rad=2.0'+ \
                       ' tmin=2020-01-01T00:01:10 tmax=2020-01-01T00:03:40'+ \
                       ' emin=0.2 emax=80.0'+ \
                       ' logfile="ctselect_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check ctselect --help
        self._check_help(ctselect)

        # Return
        return

    # Test ctselect from Python
    def _test_python(self):
        """
        Test ctselect from Python
        """
        # Allocate empty ctselect tool
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
        select['inobs']   = self._events
        select['rad']     = 2.0
        select['tmin']    = '2020-01-01T00:01:10'
        select['tmax']    = '2020-01-01T00:03:40'
        select['emin']    = 0.2
        select['emax']    = 80.0
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
        # level than before and without any selection. Since the tools
        # still holds the same selected observation container from before
        # the number of events will be identical.
        cpy_select['rad']     = 'NONE'
        cpy_select['tmin']    = 'NONE'
        cpy_select['tmax']    = 'NONE'
        cpy_select['emin']    = 'NONE'
        cpy_select['emax']    = 'NONE'
        cpy_select['outobs']  = 'ctselect_py2.fits'
        cpy_select['logfile'] = 'ctselect_py2.log'
        cpy_select['chatter'] = 3
        cpy_select['publish'] = True
        cpy_select.logFileOpen()  # Needed to get a new log file
        cpy_select.execute()

        # Check result file
        self._check_result_file('ctselect_py2.fits')

        # Execute again the copy of ctselect tool again, now manually specifying the
        # RoI centres. We also set the minimum values
        # to valid event selections, but since the maximum values are still
        # 'NONE', no event selections should occur.
        cpy_select['usepnt']   = False
        cpy_select['ra']       = 83.63
        cpy_select['dec']      = 22.01
        cpy_select['tmin']     = '2020-01-01T00:01:10'
        cpy_select['emin']     = 0.2
        cpy_select['usethres'] = 'USER'
        cpy_select['outobs']   = 'ctselect_py3.fits'
        cpy_select['logfile']  = 'ctselect_py3.log'
        cpy_select['chatter']  = 4
        cpy_select.logFileOpen()  # Needed to get a new log file
        cpy_select.execute()

        # Check result file
        self._check_result_file('ctselect_py3.fits')

        # Now clear copy of ctselect tool
        cpy_select.clear()

        # Check that cleared ctselect tool holds no observations and events
        self._check_obs(cpy_select.obs(), nobs=0)

        # Get mixed observation container
        obs = self._obs_mixed()

        # Attach response function to first observation which is the
        # event list. This is necessary to run the ctselect tool with
        # "DEFAULT" thresholds.
        obs[0].response('South_0.5h', gammalib.GCaldb('cta', 'prod2'))

        # Setup ctselect tool from observation container. An energy range
        # beyond the energies covered in the event file is specified, hence
        # an empty event list will be saved.
        select = ctools.ctselect(obs)
        select['rad']      = 2.0
        select['tmin']     = '2020-01-01T00:01:10'
        select['tmax']     = '2020-01-01T00:03:40'
        select['emin']     = 120.0
        select['emax']     = 130.0
        select['expr']     = 'DETX == 0'
        select['usethres'] = 'DEFAULT' # Has no impact as IRF has no keywords
        select['outobs']   = 'ctselect_py4.fits'
        select['logfile']  = 'ctselect_py4.log'
        select['chatter']  = 3

        # Execute tool
        select.logFileOpen()  # Needed to get a new log file
        select.execute()

        # Check result file
        self._check_result_file('ctselect_py4.fits', nevents=0)

        # Setup ctselect tool for an invalid event file
        select = ctools.ctselect()
        select['inobs']   = self._invalid_events
        select['rad']     = 2.0
        select['tmin']    = '2020-01-01T00:01:10'
        select['tmax']    = '2020-01-01T00:03:40'
        select['emin']    = 0.2
        select['emax']    = 80.0
        select['outobs']  = 'ctselect_py5.fits'
        select['logfile'] = 'ctselect_py5.log'
        select['chatter'] = 3

        # Execute tool
        select.logFileOpen()  # Needed to get a new log file
        self.test_try('Test invalid event file')
        try:
            select.execute()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Setup ctselect tool from event list with extension name. The "emax"
        # value should be ignored.
        select = ctools.ctselect()
        select['inobs']   = self._events+'[EVENTS]'
        select['rad']     = 2.0
        select['tmin']    = '2020-01-01T00:01:10'
        select['tmax']    = '2020-01-01T00:03:40'
        select['emin']    = 0.2
        select['emax']    = 0.0     # Signals that "emax" should be ignored
        select['outobs']  = 'ctselect_py6.fits'
        select['logfile'] = 'ctselect_py6.log'
        select['chatter'] = 3

        # Execute tool
        select.logFileOpen()  # Needed to get a new log file
        select.execute()

        # Check result file
        self._check_result_file('ctselect_py6.fits')

        # Now ignore the "emin" value.
        select = ctools.ctselect()
        select['inobs']   = self._events+'[EVENTS]'
        select['rad']     = 2.0
        select['tmin']    = '2020-01-01T00:01:10'
        select['tmax']    = '2020-01-01T00:03:40'
        select['emin']    = 0.0     # Signals that "emin" should be ignored
        select['emax']    = 80.0
        select['outobs']  = 'ctselect_py7.fits'
        select['logfile'] = 'ctselect_py7.log'
        select['chatter'] = 3

        # Execute tool
        select.logFileOpen()  # Needed to get a new log file
        select.execute()

        # Check result file
        self._check_result_file('ctselect_py7.fits', nevents=1504)

        # Now set "emin > emax"
        select['emin']    = 150.0
        select['emax']    = 80.0
        select['outobs']  = 'ctselect_py8.fits'
        select['logfile'] = 'ctselect_py8.log'
        select['chatter'] = 3

        # Execute tool
        select.logFileOpen()  # Needed to get a new log file
        select.execute()

        # Check result file
        self._check_result_file('ctselect_py8.fits', nevents=0)

        # Setup ctselect tool with an RoI that is displaced and of the same
        # size as the original RoI. This should cause an exception to occur
        select = ctools.ctselect()
        select['inobs']   = self._events
        select['usepnt']  = False
        select['ra']      = 83.63
        select['dec']     = 24.01
        select['rad']     = 3.0
        select['tmin']    = '2020-01-01T00:01:10'
        select['tmax']    = '2020-01-01T00:03:40'
        select['emin']    = 0.2
        select['emax']    = 80.0
        select['outobs']  = 'ctselect_py9.fits'
        select['logfile'] = 'ctselect_py9.log'
        select['chatter'] = 3

        # Execute tool
        select.logFileOpen()  # Needed to get a new log file
        self.test_try('Testing selection of invalid RoI')
        try:
            select.execute()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Now force selection and verify that the execution succeeds
        select['forcesel'] = True
        select['outobs']  = 'ctselect_py10.fits'
        select['logfile'] = 'ctselect_py10.log'

        # Execute tool
        select.logFileOpen()  # Needed to get a new log file
        select.execute()

        # Check result file
        self._check_result_file('ctselect_py10.fits', nevents=833, rad=3.)

        # Finally test a custom defined RoI
        # enclosed in the original one
        select = ctools.ctselect()
        select['inobs'] = self._events
        select['usepnt'] = False
        select['ra'] = 83.63
        select['dec'] = 24.01
        select['rad'] = 1.0
        select['tmin'] = '2020-01-01T00:01:10'
        select['tmax'] = '2020-01-01T00:03:40'
        select['emin'] = 0.2
        select['emax'] = 80.0
        select['outobs'] = 'ctselect_py11.fits'
        select['logfile'] = 'ctselect_py11.log'
        select['chatter'] = 3

        # Execute tool
        select.logFileOpen()  # Needed to get a new log file
        select.execute()

        # Check result file
        self._check_result_file('ctselect_py11.fits', nevents=117, dec=24.01, rad=1.)

        # Now test phase cut
        select = ctools.ctselect()
        select['inobs']   = self._phased_events
        select['rad']     = 2.0
        select['tmin']    = 'INDEF'
        select['tmax']    = 'INDEF'
        select['emin']    = 0.1
        select['emax']    = 100.0
        select['phase']   = '0.1:0.4,0.6:0.8'
        select['outobs']  = 'ctselect_py12.fits'
        select['logfile'] = 'ctselect_py12.log'
        select['chatter'] = 2

        # Run ctselect tool
        select.logFileOpen()   # Make sure we get a log file
        select.execute()

        # Check result file
        self._check_result_file('ctselect_py12.fits', nevents=1909, dec = 22.01)

        # Now test another phase cut
        select = ctools.ctselect()
        select['inobs']   = self._phased_events
        select['rad']     = 2.0
        select['tmin']    = 'INDEF'
        select['tmax']    = 'INDEF'
        select['emin']    = 0.1
        select['emax']    = 100.0
        select['phase']   = '0.8:0.3'
        select['outobs']  = 'ctselect_py13.fits'
        select['logfile'] = 'ctselect_py13.log'
        select['chatter'] = 2

        # Run ctselect tool
        select.logFileOpen()   # Make sure we get a log file
        select.execute()

        # Check result file
        self._check_result_file('ctselect_py13.fits', nevents=1995, dec = 22.01)

        # Test invalid minimum phase value
        self.test_try('Test invalid minimum phase value')
        try:
            select['phase']   = '-0.1:0.3'
            select['outobs']  = 'ctselect_py14.fits'
            select['logfile'] = 'ctselect_py14.log'
            select.logFileOpen()
            select.execute()
            self.test_try_failure('Exception not thrown for "-0.1:0.3"')
        except ValueError:
            self.test_try_success()

        # Test invalid minimum phase value
        self.test_try('Test invalid minimum phase value')
        try:
            select['phase']   = '1.1:0.3'
            select['outobs']  = 'ctselect_py15.fits'
            select['logfile'] = 'ctselect_py15.log'
            select.logFileOpen()
            select.execute()
            self.test_try_failure('Exception not thrown for "1.1:0.3"')
        except ValueError:
            self.test_try_success()

        # Test invalid maximum phase value
        self.test_try('Test invalid maximum phase value')
        try:
            select['phase']   = '0.1:-0.3'
            select['outobs']  = 'ctselect_py16.fits'
            select['logfile'] = 'ctselect_py16.log'
            select.logFileOpen()
            select.execute()
            self.test_try_failure('Exception not thrown for "0.1:-0.3"')
        except ValueError:
            self.test_try_success()

        # Test invalid maximum phase value
        self.test_try('Test invalid maximum phase value')
        try:
            select['phase']   = '0.1:1.1'
            select['outobs']  = 'ctselect_py17.fits'
            select['logfile'] = 'ctselect_py17.log'
            select.logFileOpen()
            select.execute()
            self.test_try_failure('Exception not thrown for "0.1:1.1"')
        except ValueError:
            self.test_try_success()

        # Test invalid phase string
        self.test_try('Test invalid phase string')
        try:
            select['phase']   = '0.1:'
            select['outobs']  = 'ctselect_py18.fits'
            select['logfile'] = 'ctselect_py18.log'
            select.logFileOpen()
            select.execute()
            self.test_try_failure('Exception not thrown for "0.1:"')
        except ValueError:
            self.test_try_success()

        # Test invalid phase string
        self.test_try('Test invalid phase string')
        try:
            select['phase']   = ':1.0'
            select['outobs']  = 'ctselect_py19.fits'
            select['logfile'] = 'ctselect_py19.log'
            select.logFileOpen()
            select.execute()
            self.test_try_failure('Exception not thrown for ":1.0"')
        except ValueError:
            self.test_try_success()

        # Test invalid phase string
        self.test_try('Test invalid phase string')
        try:
            select['phase']   = '0.1-1.0'
            select['outobs']  = 'ctselect_py20.fits'
            select['logfile'] = 'ctselect_py20.log'
            select.logFileOpen()
            select.execute()
            self.test_try_failure('Exception not thrown for "0.1-1.0"')
        except ValueError:
            self.test_try_success()

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nevents=744, dec=22.51, rad=2.0):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Event list file name
        nevents : int, optional
            Expected number of events
        dec : float, optional
            Expected Declination (deg)
        rad : float, optional
            Expected radius (deg)
        """
        # Open result file
        events = gammalib.GCTAEventList(filename)

        # Check event list
        self._check_events(events, nevents=nevents, dec=dec, rad=rad)

        # Return
        return

    # Check observation and event list
    def _check_obs(self, obs, nobs=1, nevents=744):
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
    def _check_events(self, events, nevents=744, dec=22.51, rad=2.0):
        """
        Check event list

        Parameters
        ----------
        events : `~gammalib.GCTAEventList`
            Event list
        nevents : int, optional
            Expected number of events
        dec : float, optional
            Expected Declination (deg)
        rad : float, optional
            Expected radius (deg)
        """
        # Check event list
        self.test_value(events.size(), nevents, 'Check number of events')
        self.test_value(events.roi().centre().dir().ra_deg(), 83.63,
                        'Check for RoI centre Right Ascension')
        self.test_value(events.roi().centre().dir().dec_deg(), dec,
                        'Check for RoI centre Declination')
        self.test_value(events.roi().radius(), rad,
                        'Check for RoI radius')

        # Return
        return
