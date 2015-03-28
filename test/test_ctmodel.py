#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctmodel tool.
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
        model["incube"].filename("NONE")
        model["outcube"].filename("modmap.fits")
        model["inmodel"].filename(self.model_name)
        model["inobs"].filename("NONE")
        model["expcube"].filename("NONE")
        model["psfcube"].filename("NONE")
        model["bkgcube"].filename("NONE")
        model["caldb"].string(self.caldb)
        model["irf"].string(self.irf)
        model["rad"].real(5.0)
        model["ra"].real(83.63)
        model["dec"].real(22.01)
        model["tmin"].real(0.0)
        model["tmax"].real(1800.0)
        model["emin"].real(0.1)
        model["emax"].real(100.0)
        model["enumbins"].integer(20)
        model["nxpix"].integer(200)
        model["nypix"].integer(200)
        model["binsz"].real(0.02)
        model["coordsys"].string("CEL")
        model["proj"].string("CAR")
        model["xref"].real(83.63)
        model["yref"].real(22.01)
        
        # Run tool
        self.test_try("Run ctmodel")
        try:
            model.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctmodel.")

        # Save counts cube
        self.test_try("Save model cube")
        try:
            model.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving model cube.")

        # Return
        return
