#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctobssim tool.
#
# Copyright (C) 2014 Juergen Knoedlseder
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
        self.model_name  = "data/crab.xml"
        self.events_name = "events.fits"
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
        self.name("ctobssim")

        # Append tests
        self.append(self.test_functional, "Test ctobssim functionality")

        # Return
        return

    # Test ctobssim functionnality
    def test_functional(self):
        """
        Test ctobssim functionnality.
        """
        # Set-up ctobssim
        sim = ctools.ctobssim()
        sim["infile"].filename(self.model_name)
        sim["outfile"].filename("events.fits")
        sim["caldb"].string(self.caldb)
        sim["irf"].string(self.irf)
        sim["ra"].real(83.63)
        sim["dec"].real(22.01)
        sim["rad"].real(10.0)
        sim["tmin"].real(0.0)
        sim["tmax"].real(1800.0)
        sim["emin"].real(0.1)
        sim["emax"].real(100.0)
        
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
        self.test_value(evt.size(), 4134, "Number of events is not 4134")
        
        # Save events
        self.test_try("Save events")
        try:
            sim.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving events.")
