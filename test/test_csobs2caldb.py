#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csobs2caldb script.
#
# Copyright (C) 2016 Juergen Knoedlseder
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
import cscripts


# ================================= #
# Test class for csobs2caldb script #
# ================================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for csobs2caldb script.

    This test class makes unit tests for the csobs2caldb script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set members
        self._inobs   = "data/obs_unbinned.xml"

        # Clean calibration database
        os.system("rm -rf caldb >/dev/null 2>&1")
        os.system("mkdir -p caldb >/dev/null 2>&1")

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("csobs2caldb")

        # Append tests
        self.append(self._test_cmd, "Test csobs2caldb on command line")
        self.append(self._test_python, "Test csobs2caldb from Python")

        # Return
        return

    # Test csobs2caldb on command line
    def _test_cmd(self):
        """
        Test csobs2caldb on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README"):
            csobs2caldb = "../cscripts/csobs2caldb.py"
        else:
            csobs2caldb = "csobs2caldb"

        # Setup csobs2caldb command
        cmd = csobs2caldb+' inobs="'+self._inobs+'"'+ \
                          ' caldb="cta" irf="Cmd1" rootdir="caldb"'+ \
                          ' logfile="csobs2caldb_cmd1.log" chatter=1'

        # Execute csobs2caldb, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         "Successful csobs2caldb execution on command line")

        # Check response
        self._check_response("cta", "Cmd1")

        # Setup csobs2caldb command
        cmd = csobs2caldb+' inobs="xml_file_that_does_not_exist.xml"'+ \
                          ' caldb="cta" irf="Cmd2" rootdir="caldb"'+ \
                          ' logfile="csobs2caldb_cmd2.log"'

        # Execute csobs2caldb, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         "Failure of csobs2caldb execution on command line")

        # Return
        return

    # Test csobs2caldb from Python
    def _test_python(self):
        """
        Test csobs2caldb from Python.
        """
        # Set-up csobs2caldb using "inobs" parameter
        irf = cscripts.csobs2caldb()
        irf["inobs"]   = self._inobs
        irf["caldb"]   = "cta"
        irf["irf"]     = "Py1"
        irf["rootdir"] = "caldb"
        irf["logfile"] = "csobs2caldb_py1.log"
        irf["chatter"] = 2

        # Run csobs2caldb script and save IRF
        irf.logFileOpen()   # Make sure we get a log file
        irf.run()
        irf.save()

        # Check response
        self._check_response("cta", "Py1")

        # Load observation container
        obs = gammalib.GObservations(self._inobs)

        # Set-up csobs2caldb using observation container
        irf = cscripts.csobs2caldb(obs)
        irf["caldb"]   = "cta"
        irf["irf"]     = "Py2"
        irf["rootdir"] = "caldb"
        irf["logfile"] = "csobs2caldb_py2.log"
        irf["chatter"] = 3

        # Execute csobs2caldb script
        irf.execute()

        # Check response
        self._check_response("cta", "Py2")

        # Return
        return

    # Check response
    def _check_response(self, caldb, irf):
        """
        Check response.
        """
        # Open calibration database
        db = gammalib.GCaldb("caldb")
        db.open("cta", caldb)

        # Get filenames of response components
        expr  = "NAME("+irf+")"
        aeff  = gammalib.GFilename(db.filename("","","EFF_AREA","","",expr))
        psf   = gammalib.GFilename(db.filename("","","RPSF","","",expr))
        edisp = gammalib.GFilename(db.filename("","","EDISP","","",expr))
        bgd   = gammalib.GFilename(db.filename("","","BGD","","",expr))

        # Check whether files exist
        self.test_assert(aeff.exists(), "Effective area file exists")
        self.test_assert(psf.exists(), "Point spread function file exists")
        self.test_assert(edisp.exists(), "Energy dispersion file exists")
        self.test_assert(bgd.exists(), "Background file exists")
        
        # Return
        return
