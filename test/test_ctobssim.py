#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctobssim tool
#
# Copyright (C) 2014-2017 Juergen Knoedlseder
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
                       ' outevents="ctobssim_cmd1.fits"'+ \
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
        evt = gammalib.GCTAEventList('ctobssim_cmd1.fits')
        self._test_list(evt, 7105)

        # Setup ctobssim command
        cmd = ctobssim+' inmodel="model_that_does_not_exist.xml"'+ \
                       ' outevents="ctobssim_cmd2.fits"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'" '+ \
                       ' ra=83.63 dec=22.01 rad=10.0'+ \
                       ' tmin=0.0 tmax=1800.0 emin=0.1 emax=100.0'+ \
                       ' logfile="ctobssim_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctobssim --help option
        cmd = ctobssim+' --help'

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

    # Test ctobssim from Python
    def _test_ctobssim_python(self):
        """
        Test ctobssim from Python
        """
        # Allocate ctobssim
        sim = ctools.ctobssim()

        # Check that empty ctobssim tool holds an empty observation
        self.test_value(sim.obs().size(), 0, 'Check number of observations')

        # Check that saving saves an empty model definition file
        sim['outevents'] = 'ctobssim_py0.fits'
        sim['logfile']   = 'ctobssim_py0.log'
        sim.logFileOpen()
        sim.save()
        self.test_assert(not os.path.isfile('ctobssim_py0.fits'),
             'Check that no event list has been created')

        # Check that clearing does not lead to an exception or segfault
        sim.clear()

        # Now set ctobssim parameters
        sim = ctools.ctobssim()
        sim['inmodel']   = self._model
        sim['caldb']     = self._caldb
        sim['irf']       = self._irf
        sim['ra']        = 83.63
        sim['dec']       = 22.01
        sim['rad']       = 10.0
        sim['tmin']      = 0.0
        sim['tmax']      = 1800.0
        sim['emin']      = 0.1
        sim['emax']      = 100.0
        sim['outevents'] = 'ctobssim_py1.fits'
        sim['logfile']   = 'ctobssim_py1.log'
        sim['chatter']   = 2

        # Run tool
        sim.logFileOpen()
        sim.run()

        # Check content of observation
        self._test_observation(sim)
        self._test_list(sim.obs()[0].events(), 7105)

        # Save events
        sim.save()

        # Load counts cube and check content.
        evt = gammalib.GCTAEventList('ctobssim_py1.fits')
        self._test_list(evt, 7105)

        # Set-up observation container
        pnts = [{'ra': 83.63, 'dec': 21.01},
                {'ra': 84.63, 'dec': 22.01},
                {'ra': 83.63, 'dec': 23.01},
                {'ra': 82.63, 'dec': 22.01}]
        obs = obsutils.set_obs_list(pnts, caldb=self._caldb, irf=self._irf)

        # Set-up ctobssim tool from observation container
        sim = ctools.ctobssim(obs)
        sim['inmodel']   = self._model
        sim['outevents'] = 'ctobssim_py2.xml'
        sim['logfile']   = 'ctobssim_py2.log'
        sim['chatter']   = 3

        # Double maximum event range
        max_rate = sim.max_rate()
        sim.max_rate(2*max_rate)
        self.test_value(sim.max_rate(), 2*max_rate, 1.0e-3,
                        'Check setting of maximum event rate')

        # Run ctobssim tool
        sim.logFileOpen()
        sim.run()

        # Retrieve observation and check content
        self._test_observation(sim, nobs=4, pnts=pnts)
        self._test_list(sim.obs()[0].events(), 6199)
        self._test_list(sim.obs()[1].events(), 6305)
        self._test_list(sim.obs()[2].events(), 6092)
        self._test_list(sim.obs()[3].events(), 6192)

        # Save events
        sim.save()

        # Load events
        obs = gammalib.GObservations('ctobssim_py2.xml')

        # Retrieve observation and check content
        self._test_list(obs[0].events(), 6199)
        self._test_list(obs[1].events(), 6305)
        self._test_list(obs[2].events(), 6092)
        self._test_list(obs[3].events(), 6192)

        # Copy ctobssim tool
        cpy_sim = sim.copy()

        # Retrieve observation and check content of copy
        self._test_observation(cpy_sim, nobs=4, pnts=pnts)
        self._test_list(cpy_sim.obs()[0].events(), 6199)
        self._test_list(cpy_sim.obs()[1].events(), 6305)
        self._test_list(cpy_sim.obs()[2].events(), 6092)
        self._test_list(cpy_sim.obs()[3].events(), 6192)

        # Execute copy of ctobssim tool again, now with a higher chatter
        # level than before
        cpy_sim['outevents'] = 'ctobssim_py3.xml'
        cpy_sim['logfile']   = 'ctobssim_py3.log'
        cpy_sim['chatter']   = 4
        cpy_sim['publish']   = True
        cpy_sim.logFileOpen()  # Needed to get a new log file
        cpy_sim.execute()

        # Check result file
        obs = gammalib.GObservations('ctobssim_py3.xml')
        self.test_value(obs.size(), 4, 'Check for number of observations')
        self._test_list(obs[0].events(), 6199)
        self._test_list(obs[1].events(), 6305)
        self._test_list(obs[2].events(), 6092)
        self._test_list(obs[3].events(), 6192)

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
            self.test_value(obs.livetime(), 1764.0, 1.0e-6,
                 'Livetime is 1764 sec')
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
