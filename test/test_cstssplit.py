#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cstssplit script.
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
    Test class for cstssplit script.

    This test class makes unit tests for the csstsplit script by using it
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
        self.name('cstssplit')

        # Append tests
        self.append(self._test_cmd, "Test cstssplit on command line")
        #self.append(self._test_python, "Test cstssplit from Python")

        # Return
        return

    # Test cstssplit on command line
    def _test_cmd(self):
        """
        Test cstssplit on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README.md"):
            cstssplit = "../cscripts/cstssplit.py"
        else:
            cstssplit = "cstssplit"

        # Setup cstssplit command
        cmd = cstssplit+' inobs="selected_events.fits"'+ \
                        ' inmodel="data/crab.xml"'+ \
                        ' srcname="Crab"'+ \
                        ' caldb="irf"'+ \
                        ' irf="cta_dummy_irf"'+ \
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

        # Setup cstssplit command
        cmd = cstssplit+' inobs="selected_events.fits"'+ \
                        ' inmodel="data/crab.xml"'+ \
                        ' srcname="Does_not_exist"'+ \
                        ' caldb="irf"'+ \
                        ' irf="cta_dummy_irf"'+ \
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

    # Test cstssplit from Python
    def _test_python(self):
        """
        Test cstssplit from Python.
        """
        # Set-up cstssplit
        tssplit = cscripts.cstssplit()
        tssplit['inobs'] = 'obs_selected.xml'
        tssplit['inmodel'] = 'data/crab.xml'
        tssplit['xref']   = 83.6331
        tssplit['yref']  = 22.01
        tssplit['coordsys']   = 'CEL'
        tssplit['proj'] = 'CAR'
        tssplit['binsz']  = 0.05
        tssplit['nxpix']  = 5
        tssplit['nypix']  = 5
        tssplit['compute_null']  = False
        tssplit['bins_per_job']  = 5
        tssplit['outfile']  = 'tsmap_cmd_py1.txt'


        # Run cstssplit script and save run list
        tssplit.logFileOpen()   # Make sure we get a log file
        tssplit.run()
        tssplit.save()

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
        self.test_assert(fname.exists(),
                          "Output of cstssplit exists")   
        
        # Return
        return
