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


# ============================ #
# Test class for ctobssim tool #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctobssim tool.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set members
        self.model_name = "data/crab.xml"
        self.caldb      = "irf"
        self.irf        = "cta_dummy_irf"

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
        self.append(self.test_functional, "Test ctobssim functionality")
        self.append(self.test_container, "Test ctobssim on observation container")

        # Return
        return

    # Test ctobssim functionnality
    def test_functional(self):
        """
        Test ctobssim functionnality.
        """
        # Set-up ctobssim
        sim = ctools.ctobssim()
        sim["inmodel"]   = self.model_name
        sim["outevents"] = "events.fits"
        sim["caldb"]     = self.caldb
        sim["irf"]       = self.irf
        sim["ra"]        = 83.63
        sim["dec"]       = 22.01
        sim["rad"]       = 10.0
        sim["tmin"]      = 0.0
        sim["tmax"]      = 1800.0
        sim["emin"]      = 0.1
        sim["emax"]      = 100.0

        # Run tool
        self.test_try("Run ctobssim")
        try:
            sim.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctobssim.")

        # Retrieve observation and check content
        obs = gammalib.GCTAObservation(sim.obs()[0])
        pnt = obs.pointing()
        evt = obs.events()
        self.test_assert(obs.instrument() == "CTA", "Observation not a CTA observation")
        self.test_value(sim.obs().size(), 1, "There is not a single observation")
        self.test_value(obs.ontime(), 1800.0, 1.0e-6, "Ontime is not 1800 sec")
        self.test_value(obs.livetime(), 1710.0, 1.0e-6, "Livetime is not 1710 sec")
        self.test_value(pnt.dir().ra_deg(), 83.63, 1.0e-6, "ROI Right Ascension is not 83.63 deg")
        self.test_value(pnt.dir().dec_deg(), 22.01, 1.0e-6, "ROI Declination is not 22.01 deg")
        self.test_value(evt.size(), 4119, "Number of events is not 4119")

        # Save events
        self.test_try("Save events")
        try:
            sim.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving events.")

    # Test ctobssim on observation container
    def test_container(self):
        """
        Test ctobssim on observation container.
        """
        # Set-up observation container
        obs = self.set_obs(4)

        # Set-up ctobssim
        sim = ctools.ctobssim(obs)
        sim["outevents"] = "sim_events.xml"
        sim["inmodel"]   = self.model_name
        # Run tool
        self.test_try("Run ctobssim")
        try:
            sim.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctobssim.")

        # Retrieve observation and check content
        self.test_value(sim.obs().size(), 4, "There are not 4 observations")

        # Save events
        self.test_try("Save events")
        try:
            sim.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving events.")

    # Setup observation container
    def set_obs(self, number):
        """
        Setup observation container.
        """
        # Initialise empty observation container
        obs = gammalib.GObservations()

        # Initialise first time and identifier
        tstart = 0.0
        offset = -float(number)/2.0

        # Loop over number of observations
        for i in range(number):
            id      = "%6.6d" % i
            obs_cta = self.set_one_obs(id, offset)
            offset += 1.0
            tstart += 1830.0
            obs.append(obs_cta)

        # Return observation container
        return obs

    # Setup one observation
    def set_one_obs(self, id, offset, \
                    tstart=0.0, duration=1800.0, deadc=0.95, \
                    emin=0.1, emax=100.0, rad=5.0, \
                    irf="cta_dummy_irf", caldb="dummy"):
        """
        Setup one observation for test purposes.
        """
        # Allocate CTA observation
        obs_cta = gammalib.GCTAObservation()

        # Set calibration database
        db = gammalib.GCaldb()
        if (gammalib.dir_exists(self.caldb)):
            db.rootdir(self.caldb)
        else:
            db.open("cta", self.caldb)

        # Set pointing direction
        pntdir = gammalib.GSkyDir()
        pntdir.radec_deg(83.63, 22.01+offset)
        pnt = gammalib.GCTAPointing()
        pnt.dir(pntdir)
        obs_cta.pointing(pnt)

        # Set ROI
        roi     = gammalib.GCTARoi()
        instdir = gammalib.GCTAInstDir()
        instdir.dir(pntdir)
        roi.centre(instdir)
        roi.radius(rad)

        # Set GTI
        gti = gammalib.GGti()
        gti.append(gammalib.GTime(tstart), gammalib.GTime(tstart+duration))

        # Set energy boundaries
        ebounds = gammalib.GEbounds(gammalib.GEnergy(emin, "TeV"), \
                                    gammalib.GEnergy(emax, "TeV"))

        # Allocate event list
        events = gammalib.GCTAEventList()
        events.roi(roi)
        events.gti(gti)
        events.ebounds(ebounds)
        obs_cta.events(events)

        # Set instrument response
        obs_cta.response(self.irf, db)

        # Set ontime, livetime, and deadtime correction factor
        obs_cta.ontime(duration)
        obs_cta.livetime(duration*deadc)
        obs_cta.deadc(deadc)
        obs_cta.id(id)

        # Return CTA observation
        return obs_cta
