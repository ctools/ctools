#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cspull script.
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


# ============================ #
# Test class for cspull script #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for cspull script.

    This test class makes unit tests for the cspull script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("cspull")

        # Append tests
        self.append(self._test_cmd, "Test cspull on command line")
        self.append(self._test_python, "Test cspull from Python")

        # Return
        return

    # Test cspull on command line
    def _test_cmd(self):
        """
        Test cspull on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README.md"):
            cspull = "../cscripts/cspull.py"
        else:
            cspull = "cspull"

        # Setup cspull command
        cmd = cspull+' inmodel="data/crab.xml"'+ \
                     ' outfile="cspull_cmd1.dat"'+ \
                     ' ntrials=3 caldb="prod2" irf="South_0.5h"'+ \
                     ' ra=83.6331 dec=22.0145 emin=0.1 emax=100.0'+ \
                     ' enumbins=0 tmax=1800.0 deadc=0.95 rad=5.0'+ \
                     ' npix=200 binsz=0.05'+ \
                     ' logfile="cspull_cmd1.log" chatter=1'

        # Execute cspull, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         "Successful cspull execution on command line")

        # Check pull distribution file
        self._check_pull_file("cspull_cmd1.dat")

        # Setup cspull command
        cmd = cspull+' inmodel="model_that_does_not_exist.xml"'+ \
                     ' outfile="cspull_cmd1.dat"'+ \
                     ' ntrials=3 caldb="prod2" irf="South_0.5h"'+ \
                     ' ra=83.6331 dec=22.0145 emin=0.1 emax=100.0'+ \
                     ' enumbins=0 tmax=1800.0 deadc=0.95 rad=5.0'+ \
                     ' npix=200 binsz=0.05'+ \
                     ' logfile="cspull_cmd2.log"'

        # Execute cspull, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         "Failure of cspull execution on command line")

        # Return
        return

    # Test cspull from Python
    def _test_python(self):
        """
        Test cspull from Python.
        """
        # Set-up unbinned cspull
        pull = cscripts.cspull()
        pull["inmodel"]  = "data/crab.xml"
        pull["outfile"]  = "cspull_py1.dat"
        pull["ntrials"]  = 3
        pull["caldb"]    = "prod2"
        pull["irf"]      = "South_0.5h"
        pull["ra"]       = 83.6331
        pull["dec"]      = 22.0145
        pull["emin"]     = 0.1
        pull["emax"]     = 100.0
        pull["enumbins"] = 0
        pull["tmax"]     = 1800.0
        pull["deadc"]    = 0.95
        pull["rad"]      = 5.0
        pull["logfile"]  = "cspull_py1.log"
        pull["chatter"]  = 2

        # Run cspull script
        pull.logFileOpen()   # Make sure we get a log file
        pull.run()
        #pull.save()

        # Check pull distribution file
        self._check_pull_file("cspull_py1.dat")

        # Set-up binned cspull
        pull = cscripts.cspull()
        pull["inmodel"]  = "data/crab.xml"
        pull["outfile"]  = "cspull_py2.dat"
        pull["ntrials"]  = 3
        pull["caldb"]    = "prod2"
        pull["irf"]      = "South_0.5h"
        pull["ra"]       = 83.6331
        pull["dec"]      = 22.0145
        pull["emin"]     = 0.1
        pull["emax"]     = 100.0
        pull["enumbins"] = 10
        pull["tmax"]     = 1800.0
        pull["deadc"]    = 0.95
        pull["rad"]      = 5.0
        pull["npix"]     = 100
        pull["binsz"]    = 0.02
        pull["coordsys"] = "CEL"
        pull["proj"]     = "TAN"
        pull["logfile"]  = "cspull_py2.log"
        pull["chatter"]  = 3

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file("cspull_py2.dat")

        # Set-up cspull from event list
        pull = cscripts.cspull()
        pull["inobs"]    = "data/crab_events.fits"
        pull["inmodel"]  = "data/crab.xml"
        pull["outfile"]  = "cspull_py3.dat"
        pull["ntrials"]  = 3
        pull["caldb"]    = "prod2"
        pull["irf"]      = "South_0.5h"
        pull["enumbins"] = 0
        pull["logfile"]  = "cspull_py3.log"
        pull["chatter"]  = 4

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file("cspull_py3.dat")

        # Build observation container with unbinned observation
        cta = gammalib.GCTAObservation("data/crab_events.fits")
        obs = gammalib.GObservations()
        obs.append(cta)

        # Set-up cspull from observation container with unbinned observation
        pull = cscripts.cspull(obs)
        pull["inmodel"]  = "data/crab.xml"
        pull["outfile"]  = "cspull_py4.dat"
        pull["ntrials"]  = 3
        pull["caldb"]    = "prod2"
        pull["irf"]      = "South_0.5h"
        pull["enumbins"] = 0
        pull["logfile"]  = "cspull_py4.log"
        pull["chatter"]  = 4

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file("cspull_py4.dat")

        # Set-up stacked cspull
        pull = cscripts.cspull()
        pull["inobs"]    = "data/obs_stacked.xml"
        pull["inmodel"]  = "data/crab_bkgcube.xml"
        pull["outfile"]  = "cspull_py5.dat"
        pull["ntrials"]  = 3
        pull["caldb"]    = "prod2"
        pull["irf"]      = "South_0.5h"
        pull["enumbins"] = 0
        pull["logfile"]  = "cspull_py5.log"
        pull["chatter"]  = 4

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file("cspull_py5.dat")

        # Return
        return

    # Check pull file
    def _check_pull_file(self, filename):
        """
        Check pull file.
        """
        # Open pull file as CSV file
        pulls = gammalib.GCsv(filename)

        # Check dimensions
        self.test_value(pulls.nrows(), 4,
                        'Check for 4 rows in pull file')

        # Return
        return
