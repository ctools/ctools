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
import os
import gammalib
import ctools
from cscripts import obsutils


# ============================ #
# Test class for ctobssim tool #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctobssim tools.
    
    This test class makes unit tests for the ctobssim tool by using it from
    the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set members
        self._model_name = "data/crab.xml"
        self._caldb      = "irf"
        self._irf        = "cta_dummy_irf"

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("ctobssim")

        # Append tests
        self.append(self._test_ctobssim_cmd, "Test ctobssim on command line")
        self.append(self._test_ctobssim_python, "Test ctobssim from Python")

        # Return
        return

    # Test ctobssim on command line
    def _test_ctobssim_cmd(self):
        """
        Test ctobssim on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README.md"):
            ctobssim = "../src/ctobssim/ctobssim"
        else:
            ctobssim = "ctobssim"

        # Setup ctobssim command
        cmd = ctobssim+' inmodel="data/crab.xml" outevents="events.fits"'+ \
                       ' caldb="irf" irf="cta_dummy_irf" '+ \
                       ' ra=83.63 dec=22.01 rad=10.0'+ \
                       ' tmin=0.0 tmax=1800.0 emin=0.1 emax=100.0'

        # Execute ctobssim, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         "Successful ctobssim execution on command line")

        # Load counts cube and check content.
        evt = gammalib.GCTAEventList("events.fits")
        self._test_list(evt, 4119)

        # Setup ctobssim command
        cmd = ctobssim+' inmodel="model_that_does_not_exist.xml"'+ \
                       ' outevents="events.fits"'+ \
                       ' caldb="irf" irf="cta_dummy_irf" '+ \
                       ' ra=83.63 dec=22.01 rad=10.0'+ \
                       ' tmin=0.0 tmax=1800.0 emin=0.1 emax=100.0'

        # Execute ctobssim, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         "Failure of ctobssim execution on command line")

        # Return
        return

    # Test ctobssim from Python
    def _test_ctobssim_python(self):
        """
        Test ctobssim from Python.
        """
        # Set-up ctobssim
        sim = ctools.ctobssim()
        sim["inmodel"]   = self._model_name
        sim["outevents"] = "events.fits"
        sim["caldb"]     = self._caldb
        sim["irf"]       = self._irf
        sim["ra"]        = 83.63
        sim["dec"]       = 22.01
        sim["rad"]       = 10.0
        sim["tmin"]      = 0.0
        sim["tmax"]      = 1800.0
        sim["emin"]      = 0.1
        sim["emax"]      = 100.0

        # Run tool
        sim.run()

        # Check content of observation
        self._test_observation(sim)
        self._test_list(sim.obs()[0].events(), 4119)

        # Save events
        sim.save()

        # Load counts cube and check content.
        evt = gammalib.GCTAEventList("events.fits")
        self._test_list(evt, 4119)

        # Set-up observation container
        pnts = [{'ra': 83.63, 'dec': 21.01},
                {'ra': 84.63, 'dec': 22.01},
                {'ra': 83.63, 'dec': 23.01},
                {'ra': 82.63, 'dec': 22.01}]
        obs = obsutils.set_obs_list(pnts, caldb=self._caldb, irf=self._irf)

        # Set-up ctobssim from observation container
        sim = ctools.ctobssim(obs)
        sim["outevents"] = "sim_events.xml"
        sim["inmodel"]   = self._model_name

        # Run tool
        sim.run()

        # Retrieve observation and check content
        self._test_observation(sim, nobs=4, pnts=pnts)
        self._test_list(sim.obs()[0].events(), 4063)
        self._test_list(sim.obs()[1].events(), 4059)
        self._test_list(sim.obs()[2].events(), 4010)
        self._test_list(sim.obs()[3].events(), 4138)

        # Save events
        sim.save()

        # Load events
        obs = gammalib.GObservations("sim_events.xml")

        # Retrieve observation and check content
        self._test_list(obs[0].events(), 4063)
        self._test_list(obs[1].events(), 4059)
        self._test_list(obs[2].events(), 4010)
        self._test_list(obs[3].events(), 4138)

        # Set-up ctobssim with invalid event file
        sim = ctools.ctobssim()
        sim["inmodel"]   = "model_file_that_does_not_exist.xml"
        sim["outevents"] = "events.fits"
        sim["caldb"]     = self._caldb
        sim["irf"]       = self._irf
        sim["ra"]        = 83.63
        sim["dec"]       = 22.01
        sim["rad"]       = 10.0
        sim["tmin"]      = 0.0
        sim["tmax"]      = 1800.0
        sim["emin"]      = 0.1
        sim["emax"]      = 100.0

        # Run ctbin tool
        self.test_try("Run ctobssim with invalid model file")
        try:
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
        Test content of an observation.
        
        Args:
            ctobssim: ctobssim instance

        Kwargs:
            nobs:     Number of observations
            pnts:     List of pointing dictionaries
        """
        # Test observation container
        self.test_value(ctobssim.obs().size(), nobs,
             "There is one observation")
        for i in range(ctobssim.obs().size()):
            obs = gammalib.GCTAObservation(ctobssim.obs()[i])
            pnt = obs.pointing()
            self.test_assert(obs.instrument() == "CTA",
                 "Observation is CTA observation")
            self.test_value(obs.ontime(), 1800.0, 1.0e-6,
                 "Ontime is 1800 sec")
            self.test_value(obs.livetime(), 1710.0, 1.0e-6,
                 "Livetime is 1710 sec")
            self.test_value(pnt.dir().ra_deg(), pnts[i]['ra'], 1.0e-6,
                 "Pointing Right Ascension is "+str(pnts[i]['ra'])+" deg")
            self.test_value(pnt.dir().dec_deg(), pnts[i]['dec'], 1.0e-6,
                 "Pointing Declination is "+str(pnts[i]['dec'])+" deg")

        # Return
        return

    # Check event list
    def _test_list(self, list, nevents):
        """
        Test content of event list."
        
        Args:
            list:    Event list
            nevents: Expected number of events.
        """
        # Test event list
        self.test_value(list.size(), nevents, str(nevents)+" elements")
        self.test_value(list.number(), nevents, str(nevents)+" events")

        # Return
        return
