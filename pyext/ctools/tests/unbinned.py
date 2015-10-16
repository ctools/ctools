# ==========================================================================
# This script tests the unbinned analysis.
#
# Copyright (C) 2015 Christoph Deil
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
from tempfile import mkdtemp
import unittest
import ctools


# ================== #
# Analysis pipelines #
# ================== #
class pipelines(unittest.TestCase):
    """
    Test a simple unbinned analysis pipeline.
    """
    # Setup of test class
    def setUp(self):
        gammalib_dir    = os.environ['GAMMALIB']
        self.model_name = os.path.join(gammalib_dir, "share/models/crab.xml")
        self.workdir    = mkdtemp()

    # Test unbinned analysis with FITS files
    def test_unbinned_fits(self):
        """
        Test unbinned pipeline with FITS file saving.
        """
        # Set script parameters
        model_name           = self.model_name
        events_name          = os.path.join(self.workdir, "events.fits")
        selected_events_name = os.path.join(self.workdir, "selected_events.fits")
        result_name          = os.path.join(self.workdir, "results.xml")
        caldb                = "prod2"
        irf                  = "South_50h"
        ra                   =   83.63
        dec                  =   22.01
        rad_sim              =   10.0
        tstart               =    0.0
        tstop                = 1800.0
        emin                 =    0.1
        emax                 =  100.0
        rad_select           =    3.0

        # Simulate events
        sim = ctools.ctobssim()
        sim["inmodel"]   = model_name
        sim["outevents"] = events_name
        sim["caldb"]     = caldb
        sim["irf"]       = irf
        sim["ra"]        = ra
        sim["dec"]       = dec
        sim["rad"]       = rad_sim
        sim["tmin"]      = tstart
        sim["tmax"]      = tstop
        sim["emin"]      = emin
        sim["emax"]      = emax
        sim.execute()

        # Select events
        select = ctools.ctselect()
        select["inobs"]  = events_name
        select["outobs"] = selected_events_name
        select["ra"]     = ra
        select["dec"]    = dec
        select["rad"]    = rad_select
        select["tmin"]   = tstart
        select["tmax"]   = tstop
        select["emin"]   = emin
        select["emax"]   = emax
        select.execute()

        # Perform maximum likelihood fitting
        like = ctools.ctlike()
        like["inobs"]    = selected_events_name
        like["inmodel"]  = model_name
        like["outmodel"] = result_name
        like["caldb"]    = caldb
        like["irf"]      = irf
        like.execute()

        # TODO: add asserts
        assert True

    # Test unbinned analysis in memory
    def test_unbinned_mem(self):
        """
        Test unbinned in-memory pipeline.
        """
        # Set script parameters
        model_name = self.model_name
        caldb      = "prod2"
        irf        = "South_50h"
        ra         =   83.63
        dec        =   22.01
        rad_sim    =   10.0
        tstart     =    0.0
        tstop      = 1800.0
        emin       =    0.1
        emax       =  100.0
        rad_select =    3.0

        # Simulate events
        sim = ctools.ctobssim()
        sim["inmodel"] = model_name
        sim["caldb"]   = caldb
        sim["irf"]     = irf
        sim["ra"]      = ra
        sim["dec"]     = dec
        sim["rad"]     = rad_sim
        sim["tmin"]    = tstart
        sim["tmax"]    = tstop
        sim["emin"]    = emin
        sim["emax"]    = emax
        sim.run()

        # Select events
        select = ctools.ctselect(sim.obs())
        select["ra"]   = ra
        select["dec"]  = dec
        select["rad"]  = rad_select
        select["tmin"] = tstart
        select["tmax"] = tstop
        select["emin"] = emin
        select["emax"] = emax
        select.run()

        # Perform maximum likelihood fitting
        like = ctools.ctlike(select.obs())
        like.run()

        # TODO: add asserts
        assert True
