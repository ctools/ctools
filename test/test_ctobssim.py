#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctobssim tool
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
                       ' ra=83.63 dec=22.01 rad=3.0'+ \
                       ' tmin="2020-01-01T00:00:00"'+ \
                       ' tmax="2020-01-01T00:05:00"'+ \
                       ' emin=0.1 emax=100.0'+ \
                       ' logfile="ctobssim_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Load counts cube and check content.
        evt = gammalib.GCTAEventList('ctobssim_cmd1.fits')
        self._test_list(evt, 3775)

        # Setup ctobssim command
        cmd = ctobssim+' inmodel="model_that_does_not_exist.xml"'+ \
                       ' outevents="ctobssim_cmd2.fits"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'" '+ \
                       ' ra=83.63 dec=22.01 rad=3.0'+ \
                       ' tmin="2020-01-01T00:00:00"'+ \
                       ' tmax="2020-01-01T00:05:00"'+ \
                       ' emin=0.1 emax=100.0'+ \
                       ' logfile="ctobssim_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check ctobssim --help
        self._check_help(ctobssim)

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
        sim['rad']       = 3.0
        sim['tmin']      = '2020-01-01T00:00:00'
        sim['tmax']      = '2020-01-01T00:05:00'
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
        self._test_list(sim.obs()[0].events(), 3775)

        # Save events
        sim.save()

        # Load counts cube and check content.
        evt = gammalib.GCTAEventList('ctobssim_py1.fits')
        self._test_list(evt, 3775)

        # Copy ctobssim tool
        cpy_sim = sim.copy()

        # Retrieve observation and check content of copy
        self._test_observation(cpy_sim)
        self._test_list(cpy_sim.obs()[0].events(), 3775)

        # Execute copy of ctobssim tool again, now with a higher chatter
        # level than before
        cpy_sim['outevents'] = 'ctobssim_py2.fits'
        cpy_sim['logfile']   = 'ctobssim_py2.log'
        cpy_sim['chatter']   = 3
        cpy_sim['publish']   = True
        cpy_sim.logFileOpen()  # Needed to get a new log file
        cpy_sim.execute()

        # Load counts cube and check content.
        evt = gammalib.GCTAEventList('ctobssim_py2.fits')
        self._test_list(evt, 3775)

        # Set-up observation container
        pnts = [{'ra': 83.63, 'dec': 21.01},
                {'ra': 84.63, 'dec': 22.01}]
        obs = obsutils.set_obs_list(pnts, caldb=self._caldb, irf=self._irf,
                                    duration=300.0, rad=3.0)

        # Set-up ctobssim tool from observation container
        sim = ctools.ctobssim(obs)
        sim['inmodel']   = self._model
        sim['outevents'] = 'ctobssim_py3.xml'
        sim['logfile']   = 'ctobssim_py3.log'
        sim['chatter']   = 3

        # Double maximum event range
        max_rate = sim.max_rate()
        sim.max_rate(2.0*max_rate)
        self.test_value(sim.max_rate(), 2.0*max_rate, 1.0e-3,
                        'Check setting of maximum event rate')

        # Run ctobssim tool
        sim.logFileOpen()
        sim.run()

        # Retrieve observation and check content
        self._test_observation(sim, nobs=2, pnts=pnts)
        self._test_list(sim.obs()[0].events(), 3569)
        self._test_list(sim.obs()[1].events(), 3521)

        # Save events
        sim.save()

        # Load events
        obs = gammalib.GObservations('ctobssim_py3.xml')

        # Retrieve observation and check content
        self._test_list(obs[0].events(), 3569)
        self._test_list(obs[1].events(), 3521)

        # Return
        return

    # Check observation
    def _test_observation(self, obssim, nobs=1,
                          pnts=[{'ra': 83.63, 'dec': 22.01}],
                          ontime=300.0, livetime=294.0):
        """
        Test content of an observation
        
        Parameters
        ----------
        obssim : `~ctools.ctobssim`
            ctobssim instance
        nobs : int
            Number of observations
        pnts : list of dict
            List of pointing dictionaries
        """
        # Test observation container
        self.test_value(obssim.obs().size(), nobs,
                        'Check number of observations')
        for i, obs in enumerate(obssim.obs()):
            pnt = obs.pointing().dir()
            ra  = pnts[i]['ra']
            dec = pnts[i]['dec']
            self.test_assert(obs.instrument() == 'CTA',
                 'Check if observation is CTA observation')
            self.test_value(obs.ontime(), ontime, 1.0e-6,
                 'Check if ontime is %f sec' % ontime)
            self.test_value(obs.livetime(), livetime, 1.0e-6,
                 'Check if livetime is %f sec' % livetime)
            self.test_value(pnt.ra_deg(), ra, 1.0e-6,
                 'Check if pointing Right Ascension is %f deg' % ra)
            self.test_value(pnt.dec_deg(), dec, 1.0e-6,
                 'Check if pointing Declination is %f deg' % dec)

        # Return
        return

    # Check event list
    def _test_list(self, eventlist, nevents):
        """
        Test content of event list
        
        Parameters
        ----------
        eventlist : `~gammalib.GCTAEventList`
            Event list
        nevents : int
            Expected number of events
        """
        # Test event list
        self.test_value(eventlist.size(), nevents, str(nevents)+' elements')
        self.test_value(eventlist.number(), nevents, str(nevents)+' events')

        # Return
        return
