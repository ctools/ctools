#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctedispcube tool.
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


# ============================= #
# Test class for ctedispcube tool #
# ============================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctedispcube tool.
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
        self.name("ctedispcube")

        # Append tests
        self.append(self.test_functional, "Test ctedispcube functionality")

        # Return
        return

    # Test ctedispcube functionnality
    def test_functional(self):
        """
        Test ctedispcube functionnality.
        """
        # Set-up ctedispcube
        edispcube = ctools.ctedispcube()
        edispcube["inobs"]    = self.events_name
        edispcube["incube"]   = "NONE"
        edispcube["outcube"]  = "edispcube.fits"
        edispcube["caldb"]    = self.caldb
        edispcube["irf"]      = self.irf
        edispcube["ebinalg"]  = "LOG"
        edispcube["emin"]     = 0.1
        edispcube["emax"]     = 100
        edispcube["enumbins"] = 20
        edispcube["nxpix"]    = 10
        edispcube["nypix"]    = 10
        edispcube["binsz"]    = 0.4
        edispcube["coordsys"] = "CEL"
        edispcube["proj"]     = "CAR"
        edispcube["xref"]     = 83.63
        edispcube["yref"]     = 22.01
        edispcube["mmax"]     = 0.3
        edispcube["migrabins"] = 10

        # Run tool
        self.test_try("Run ctedispcube")
        try:
            edispcube.run()
            self.test_try_success()
        #except Exception as e:
        #    msg = "Exception occured in ctedispcube: %s." % (e,)
        except:
            msg = "Exception occured in ctedispcube."
            self.test_try_failure(msg)

        # Save EDISP cube
        self.test_try("Save EDISP cube")
        try:
            edispcube.save()
            self.test_try_success()
        #except Exception as e:
        #    msg = "Exception occured in saving EDISP cube: %s." % (e,)
        except:
            msg = "Exception occured in saving EDISP cubes."
            self.test_try_failure(msg)

        # Return
        return
