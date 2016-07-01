#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cstsmapsplit script.
#
# Copyright (C) 2016 Michael Mayer
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


# =============================== #
# Test class for cstssplit script #
# =============================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for cstssplit script

    This test class makes unit tests for the csstsplit script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('cstsmapsplit')

        # Append tests
        self.append(self._test_cmd, 'Test cstsmapsplit on command line')
        #self.append(self._test_python, 'Test cstsmapsplit from Python')

        # Return
        return

    # Test cstsmapsplit on command line
    def _test_cmd(self):
        """
        Test cstsmapsplit on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README.md"):
            cstsmapsplit = "../cscripts/cstsmapsplit.py"
        else:
            cstsmapsplit = "cstsmapsplit"

        # Setup cstsmapsplit command
        cmd = cstsmapsplit+' inobs="selected_events.fits"'+ \
                           ' inmodel="data/crab.xml"'+ \
                           ' srcname="Crab"'+ \
                           ' caldb="prod2"'+ \
                           ' irf="South_0.5h"'+ \
                           ' outmap="tsmap.fits"'+ \
                           ' nxpix=5'+ \
                           ' nypix=5'+ \
                           ' binsz=0.05'+ \
                           ' coordsys="CEL"'+ \
                           ' xref=83.63'+ \
                           ' yref=22.01'+ \
                           ' proj="CAR"'+ \
                           ' compute_null=yes'+ \
                           ' bins_per_job=5'+ \
                           ' outfile="tsmap_cmds1.txt"'


        # Execute cstssplit, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         "Successful cstsssplit execution on command line")

        # Check observation definition XML file
        self._check_cmdfile('tsmap_cmds1.txt')

        # Setup cstsmapsplit command
        cmd = cstsmapsplit+' inobs="selected_events.fits"'+ \
                           ' inmodel="data/crab.xml"'+ \
                           ' srcname="Does_not_exist"'+ \
                           ' caldb="prod2"'+ \
                           ' irf="South_0.5h"'+ \
                           ' outmap="tsmap.fits"'+ \
                           ' nxpix=5'+ \
                           ' nypix=5'+ \
                           ' binsz=0.05'+ \
                           ' coordsys="CEL"'+ \
                           ' xref=83.63'+ \
                           ' yref=22.01'+ \
                           ' proj="CAR"'+ \
                           ' compute_null=yes'+ \
                           ' bins_per_job=5'+ \
                           ' outfile="tsmap_cmds2.txt"'
                        

        # Execute cstssplit, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         "Failure of cstssplit execution on command line")

        # Return
        return

    # Test cstsmapsplit from Python
    def _test_python(self):
        """
        Test cstsmapsplit from Python.
        """
        # Set-up cstsmapsplit
        tsmapsplit = cscripts.cstsmapsplit()
        tsmapsplit['inobs']        = 'obs_selected.xml'
        tsmapsplit['inmodel']      = 'data/crab.xml'
        tsmapsplit['xref']         = 83.6331
        tsmapsplit['yref']         = 22.01
        tsmapsplit['coordsys']     = 'CEL'
        tsmapsplit['proj']         = 'CAR'
        tsmapsplit['binsz']        = 0.05
        tsmapsplit['nxpix']        = 5
        tsmapsplit['nypix']        = 5
        tsmapsplit['compute_null'] = False
        tsmapsplit['bins_per_job'] = 5
        tsmapsplit['outfile']      = 'tsmap_cmd_py1.txt'

        # Run cstssplit script and save run list
        tsmapsplit.logFileOpen()   # Make sure we get a log file
        tsmapsplit.run()
        tsmapsplit.save()

        # Check observation definition XML file
        self._check_cmdfile('tsmap_cmd_py1.txt')
        
        # Return
        return

    # Check observation definition XML file
    def _check_cmdfile(self, filename):
        """
        Check ascii files containing commands
        """
        
        # Create file name instance
        fname = gammalib.GFilename(filename)
        
        # Check if execution was successful
        self.test_assert(fname.exists(), 'Check of output filename exists')
        
        # Return
        return
