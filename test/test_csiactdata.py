#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csiactdata script.
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


# ================================ #
# Test class for csiactdata script #
# ================================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for csiactdata script.

    This test class makes unit tests for the csiactdata script by using it
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
        self.name("csiactdata")

        # Append tests
        self.append(self._test_cmd, "Test csiactdata on command line")
        self.append(self._test_python, "Test csiactdata from Python")

        # Return
        return

    # Test csiactdata on command line
    def _test_cmd(self):
        """
        Test csiactdata on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README.md"):
            csiactdata = "../cscripts/csiactdata.py"
        else:
            csiactdata = "csiactdata"

        # Setup csiactdata command
        cmd = csiactdata+' datapath="iactdata"'+ \
                         ' logfile="csiactdata_cmd1.log" chatter=1'

        # Execute csiactdata, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         "Successful csiactdata execution on command line")

        # Setup csiactdata command
        cmd = csiactdata+' datapath="data_path_that_does_not_exist"'+ \
                         ' logfile="csiactdata_cmd2.log"'

        # Execute csiactdata, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         "Failure of csiactdata execution on command line")

        # Return
        return

    # Test csiactdata from Python
    def _test_python(self):
        """
        Test csiactdata from Python.
        """
        # Set-up csiactdata
        iactdata = cscripts.csiactdata()
        iactdata["datapath"] = "iactdata"
        iactdata["logfile"]  = "csiactdata_py1.log"
        iactdata["chatter"]  = 2

        # Run csiactdata script and save run list
        iactdata.logFileOpen()   # Make sure we get a log file
        iactdata.run()

        # Return
        return
