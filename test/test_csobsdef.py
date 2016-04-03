#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csobsdef script.
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


# ============================== #
# Test class for csobsdef script #
# ============================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for csobsdef script.

    This test class makes unit tests for the csobsdef script by using it
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
        self.name("csobsdef")

        # Append tests
        self.append(self._test_cmd, "Test csobsdef on command line")
        self.append(self._test_python, "Test csobsdef from Python")

        # Return
        return

    # Test csobsdef on command line
    def _test_cmd(self):
        """
        Test csobsdef on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README"):
            csobsdef = "../cscripts/csobsdef.py"
        else:
            csobsdef = "csobsdef"

        # Setup csobsdef command
        cmd = csobsdef+' inpnt="data/pntdef.dat"'+ \
                       ' outobs="obsdef_cmd1.xml"'+ \
                       ' caldb="prod2" irf="South_0.5h"'+ \
                       ' emin=0.1 emax=100.0 duration=1800.0 rad=5.0'+ \
                       ' logfile="csobsdef_cmd1.log" chatter=1'

        # Execute csobsdef, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         "Successful csobsdef execution on command line")

        # Check observation definition file
        self._check_obsdef_file("obsdef_cmd1.xml")

        # Setup csobsdef command
        cmd = csobsdef+' inpnt="obs_definition_that_does_not_exist.dat"'+ \
                       ' outobs="obsdef_cmd1.xml"'+ \
                       ' caldb="prod2" irf="South_0.5h"'+ \
                       ' emin=0.1 emax=100.0 duration=1800.0 rad=5.0'+ \
                       ' logfile="csobsdef_cmd2.log"'

        # Execute csobsdef, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         "Failure of csobsdef execution on command line")

        # Return
        return

    # Test csobsdef from Python
    def _test_python(self):
        """
        Test csobsdef from Python.
        """
        # Set-up csobsdef
        obsdef = cscripts.csobsdef()
        obsdef["inpnt"]    = "data/pntdef.dat"
        obsdef["outobs"]   = "obsdef_py1.xml"
        obsdef["caldb"]    = "prod2"
        obsdef["irf"]      = "South_0.5h"
        obsdef["emin"]     = 0.1
        obsdef["emax"]     = 100.0
        obsdef["duration"] = 1800.0
        obsdef["rad"]      = 5.0
        obsdef["logfile"]  = "csobsdef_py1.log"
        obsdef["chatter"]  = 2

        # Run csobsdef script and save XML file
        obsdef.logFileOpen()   # Make sure we get a log file
        obsdef.run()
        obsdef.save()

        # Check model file
        self._check_obsdef_file("obsdef_py1.xml")

        # Set-up csobsdef
        obsdef = cscripts.csobsdef()
        obsdef["inpnt"]    = "data/pntdef.dat"
        obsdef["outobs"]   = "obsdef_py2.xml"
        obsdef["caldb"]    = "prod2"
        obsdef["irf"]      = "South_0.5h"
        obsdef["emin"]     = 0.1
        obsdef["emax"]     = 100.0
        obsdef["duration"] = 1800.0
        obsdef["rad"]      = 5.0
        obsdef["logfile"]  = "csobsdef_py2.log"
        obsdef["chatter"]  = 3

        # Execute csobsdef script
        obsdef.execute()

        # Check model file
        self._check_obsdef_file("obsdef_py2.xml")

        # Return
        return

    # Check observation definition file
    def _check_obsdef_file(self, filename):
        """
        Check observation definition file.
        """
        # Load observation definition file
        obs = gammalib.GObservations(filename)

        # Check number of observations
        self.test_value(obs.size(), 3,
             'Check for 3 observations in observation definition file')

        # Loop over all observations
        for o in obs:
            self.test_assert(o.instrument() == "CTA",
                            'Check for "CTA" instrument')
            self.test_value(o.ontime(), 1800.0, 1.0e-6,
                            'Check for ontime of 1800 sec')
            self.test_value(o.deadc(), 0.95, 1.0e-6,
                            'Check for deadtime correction of 0.95')
            self.test_value(o.events().ebounds().emin().TeV(), 0.1, 1.0e-6,
                            'Check for minimum energy of 0.1 TeV')
            self.test_value(o.events().ebounds().emax().TeV(), 100.0, 1.0e-6,
                            'Check for maximum energy of 100 TeV')
            self.test_value(o.roi().radius(), 5.0, 1.0e-6,
                            'Check for ROI radois of 5 deg')

        # Return
        return
