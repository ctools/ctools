#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctexpcube tool.
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


# ============================= #
# Test class for ctexpcube tool #
# ============================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctexpcube tool.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set members
        self.events_name = "data/crab_events.fits"
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
        self.name("ctexpcube")

        # Append tests
        self.append(self.test_functional, "Test ctexpcube functionality")

        # Return
        return

    # Test ctexpcube functionnality
    def test_functional(self):
        """
        Test ctexpcube functionnality.
        """
        # Set-up ctexpcube
        expcube = ctools.ctexpcube()
        expcube["infile"].filename(self.events_name)
        expcube["cntmap"].filename("NONE")
        expcube["outfile"].filename("expcube.fits")
        expcube["caldb"].string(self.caldb)
        expcube["irf"].string(self.irf)
        expcube["ebinalg"].string("LOG")
        expcube["emin"].real(0.1)
        expcube["emax"].real(100.0)
        expcube["enumbins"].integer(20)
        expcube["nxpix"].integer(200)
        expcube["nypix"].integer(200)
        expcube["binsz"].real(0.02)
        expcube["coordsys"].string("CEL")
        expcube["proj"].string("CAR")
        expcube["xref"].real(83.63)
        expcube["yref"].real(22.01)
        
        # Run tool
        self.test_try("Run ctexpcube")
        try:
            expcube.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in ctexpcube.")

        # Save exposure cube
        self.test_try("Save exposure cube")
        try:
            expcube.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving exposure cube.")

        # Return
        return
