#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctobssim tool.
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
from cscripts import obsutils
from testing import test


# ============================ #
# Test class for ctobssim tool #
# ============================ #
class Test(test):
    """
    Test class for ctobssim tools
    
    This test class makes unit tests for the ctobssim tool by using it from
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
        self.name('ctobssim')

        # Append tests
        self.append(self._test_ctobssim_cmd, 'Test ctobssim on command line')
        self.append(self._test_ctobssim_python, 'Test ctobssim from Python')

        # Return
        return

    # Test ctobssim on command line
    def _test_ctobssim_cmd(self):
        """
        Test ctobssim on the command line
        """
        # Set tool name
        ctobssim = self._tool('ctobssim')

        # Setup ctobssim command
        cmd = ctobssim+' inmodel="'+self._model+'" '+ \
                       ' outevents="events.fits"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'" '+ \
                       ' ra=83.63 dec=22.01 rad=10.0'+ \
                       ' tmin=0.0 tmax=1800.0 emin=0.1 emax=100.0'+ \
                       ' logfile="ctobssim_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Load counts cube and check content.
        evt = gammalib.GCTAEventList('events.fits')
        self._test_list(evt, 6881)

        # Setup ctobssim command
        cmd = ctobssim+' inmodel="model_that_does_not_exist.xml"'+ \
                       ' outevents="events.fits"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'" '+ \
                       ' ra=83.63 dec=22.01 rad=10.0'+ \
                       ' tmin=0.0 tmax=1800.0 emin=0.1 emax=100.0'+ \
                       ' logfile="ctobssim_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctobssim from Python
    def _test_ctobssim_python(self):
        """
        Test ctobssim from Python
        """
        # Set-up ctobssim
        sim = ctools.ctobssim()
        sim['inmodel']   = self._model
        sim['outevents'] = 'events.fits'
        sim['caldb']     = self._caldb
        sim['irf']       = self._irf
        sim['ra']        = 83.63
        sim['dec']       = 22.01
        sim['rad']       = 10.0
        sim['tmin']      = 0.0
        sim['tmax']      = 1800.0
        sim['emin']      = 0.1
        sim['emax']      = 100.0
        sim['logfile']   = 'ctobssim_py1.log'
        sim['chatter']   = 2

        # Run tool
        sim.logFileOpen()
        sim.run()

        # Check content of observation
        self._test_observation(sim)
        self._test_list(sim.obs()[0].events(), 6881)

        # Save events
        sim.save()

        # Load counts cube and check content.
        evt = gammalib.GCTAEventList('events.fits')
        self._test_list(evt, 6881)

        # Set-up observation container
        pnts = [{'ra': 83.63, 'dec': 21.01},
                {'ra': 84.63, 'dec': 22.01},
                {'ra': 83.63, 'dec': 23.01},
                {'ra': 82.63, 'dec': 22.01}]
        obs = obsutils.set_obs_list(pnts, caldb=self._caldb, irf=self._irf)

        # Set-up ctobssim from observation container
        sim = ctools.ctobssim(obs)
        sim['outevents'] = 'sim_events.xml'
        sim['inmodel']   = self._model
        sim['logfile']   = 'ctobssim_py2.log'
        sim['chatter']   = 3

        # Run tool
        sim.logFileOpen()
        sim.run()

        # Retrieve observation and check content
        self._test_observation(sim, nobs=4, pnts=pnts)
        self._test_list(sim.obs()[0].events(), 6003)
        self._test_list(sim.obs()[1].events(), 6084)
        self._test_list(sim.obs()[2].events(), 5955)
        self._test_list(sim.obs()[3].events(), 6030)

        # Save events
        sim.save()

        # Load events
        obs = gammalib.GObservations('sim_events.xml')

        # Retrieve observation and check content
        self._test_list(obs[0].events(), 6003)
        self._test_list(obs[1].events(), 6084)
        self._test_list(obs[2].events(), 5955)
        self._test_list(obs[3].events(), 6030)

        # Set-up ctobssim with invalid event file
        sim = ctools.ctobssim()
        sim['inmodel']   = 'model_file_that_does_not_exist.xml'
        sim['outevents'] = 'events.fits'
        sim['caldb']     = self._caldb
        sim['irf']       = self._irf
        sim['ra']        = 83.63
        sim['dec']       = 22.01
        sim['rad']       = 10.0
        sim['tmin']      = 0.0
        sim['tmax']      = 1800.0
        sim['emin']      = 0.1
        sim['emax']      = 100.0
        sim['logfile']   = 'ctobssim_py3.log'
        sim['chatter']   = 4

        # Run ctbin tool
        self.test_try('Run ctobssim with invalid model file')
        try:
            sim.logFileOpen()
            sim.run()
            self.test_try_failure()
        except:
            self.test_try_success()

        # Return
        return

    # Check observation
    def _test_observation(self, ctobssim, nobs=1,
                          pnts=[{'ra': 83.63, 'dec': 22.01}]):
        """
        Test content of an observation
        
        Parameters
        ----------
        ctobssim : `~ctools.ctobssim`
            ctobssim instance
        nobs : int
            Number of observations
        pnts : list of dict
            List of pointing dictionaries
        """
        # Test observation container
        self.test_value(ctobssim.obs().size(), nobs,
             'There is one observation')
        for i in range(ctobssim.obs().size()):
            obs = gammalib.GCTAObservation(ctobssim.obs()[i])
            pnt = obs.pointing()
            self.test_assert(obs.instrument() == 'CTA',
                 'Observation is CTA observation')
            self.test_value(obs.ontime(), 1800.0, 1.0e-6,
                 'Ontime is 1800 sec')
            self.test_value(obs.livetime(), 1710.0, 1.0e-6,
                 'Livetime is 1710 sec')
            self.test_value(pnt.dir().ra_deg(), pnts[i]['ra'], 1.0e-6,
                 'Pointing Right Ascension is '+str(pnts[i]['ra'])+' deg')
            self.test_value(pnt.dir().dec_deg(), pnts[i]['dec'], 1.0e-6,
                 'Pointing Declination is '+str(pnts[i]['dec'])+' deg')

        # Return
        return

    # Check event list
    def _test_list(self, list, nevents):
        """
        Test content of event list
        
        Parameters
        ----------
        list : `~gammalib.GCTAEventList`
            Event list
        nevents : int
            Expected number of events
        """
        # Test event list
        self.test_value(list.size(), nevents, str(nevents)+' elements')
        self.test_value(list.number(), nevents, str(nevents)+' events')

        # Return
        return
