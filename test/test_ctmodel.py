#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctmodel tool.
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


# =========================== #
# Test class for ctmodel tool #
# =========================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctmodel tool.
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
        self.name("ctmodel")

        # Append tests
        self.append(self.test_functional, "Test ctmodel functionality")

        # Return
        return

    # Test ctmodel functionnality
    def test_functional(self):
        """
        Test ctmodel functionnality.
        """
        # Set-up ctmodel from scratch
        model = ctools.ctmodel()
        model["incube"]   = "NONE"
        model["outcube"]  = "modmap.fits"
        model["inmodel"]  = self.model_name
        model["inobs"]    = "NONE"
        model["expcube"]  = "NONE"
        model["psfcube"]  = "NONE"
        model["bkgcube"]  = "NONE"
        model["caldb"]    = self.caldb
        model["irf"]      = self.irf
        model["rad"]      = 5
        model["ra"]       = 83.63
        model["dec"]      = 22.01
        model["tmin"]     = 0
        model["tmax"]     = 1800
        model["emin"]     = 0.1
        model["emax"]     = 100
        model["enumbins"] = 20
        model["nxpix"]    = 200
        model["nypix"]    = 200
        model["binsz"]    = 0.02
        model["coordsys"] = "CEL"
        model["proj"]     = "CAR"
        model["xref"]     = 83.63
        model["yref"]     = 22.01

        # Run tool
        self.test_try("Run ctmodel")
        try:
            model.run()
            self.test_try_success()
        except Exception, e:
            msg = "Exception occured in ctmodel: %s." % (e,)
            self.test_try_failure(msg)

        # Save counts cube
        self.test_try("Save model cube")
        try:
            model.save()
            self.test_try_success()
        except Exception, e:
            msg = "Exception occured in saving model cube: %s." % (e,)
            self.test_try_failure(msg)

        # Return
        return
