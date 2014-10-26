#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cttsmap tool.
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


# ========================== #
# Test class for cttsmap tool #
# ========================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for cttsmap tool.
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
        self.model_name  = "data/crab.xml"
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
        self.name("cttsmap")

        # Append tests
        self.append(self.test_functional, "Test cttsmap functionality")

        # Return
        return

    # Test cttsmap functionnality
    def test_functional(self):
        """
        Test cttsmap functionnality.
        """
        # Set-up cttsmap
        tsmap = ctools.cttsmap()
        tsmap["debug"].boolean(True)
        tsmap["infile"].filename(self.events_name)
        tsmap["srcmdl"].filename(self.model_name)
        tsmap["srcname"].string("Crab")
        tsmap["caldb"].string(self.caldb)
        tsmap["irf"].string(self.irf)
        tsmap["outfile"].filename("tsmap.fits")
        tsmap["nxpix"].integer(5)
        tsmap["nypix"].integer(5)
        tsmap["binsz"].real(0.02)
        tsmap["coordsys"].string("CEL")
        tsmap["proj"].string("CAR")
        tsmap["xref"].real(83.63)
        tsmap["yref"].real(22.01)

        # Run tool
        self.test_try("Run cttsmap")
        try:
            tsmap.run()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in cttsmap.")

        # Save results
        self.test_try("Save results")
        try:
            tsmap.save()
            self.test_try_success()
        except:
            self.test_try_failure("Exception occured in saving results.")
        
        # Return
        return
