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
import gammalib
import cscripts
from testing import test


# ================================== #
# Test class for cstsmapsplit script #
# ================================== #
class Test(test):
    """
    Test class for cstsmapsplit script

    This test class makes unit tests for the cstsmapsplit script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

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
        self.append(self._test_python, 'Test cstsmapsplit from Python')

        # Return
        return

    # Test cstsmapsplit on command line
    def _test_cmd(self):
        """
        Test cstsmapsplit on the command line
        """
        # Set script name
        cstsmapsplit = self._script('cstsmapsplit')

        # Setup cstsmapsplit command
        cmd = cstsmapsplit+' inobs="'+self._events+'"'+ \
                           ' inmodel="'+self._model+'"'+ \
                           ' srcname="Crab"'+ \
                           ' caldb="'+self._caldb+'"'+ \
                           ' irf="'+self._irf+'"'+ \
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
                           ' outfile="cstsmapsplit_cmd1.dat"'+ \
                           ' logfile="cstsmapsplit_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check command file
        self._check_cmdfile('cstsmapsplit_cmd1.dat')

        # Setup cstsmapsplit command
        cmd = cstsmapsplit+' inobs="'+self._events+'"'+ \
                           ' inmodel="'+self._model+'"'+ \
                           ' srcname="Does_not_exist"'+ \
                           ' caldb="'+self._caldb+'"'+ \
                           ' irf="'+self._irf+'"'+ \
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
                           ' outfile="cstsmapsplit_cmd2.dat"'+ \
                           ' logfile="cstsmapsplit_cmd2.log" chatter=2'
        
        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test cstsmapsplit from Python
    def _test_python(self):
        """
        Test cstsmapsplit from Python
        """
        # Set-up cstsmapsplit
        tsmapsplit = cscripts.cstsmapsplit()
        tsmapsplit['inobs']        = self._events
        tsmapsplit['inmodel']      = self._model
        tsmapsplit['srcname']      = 'Crab'
        tsmapsplit['caldb']        = self._caldb
        tsmapsplit['irf']          = self._irf
        tsmapsplit['outmap']       = 'tsmap.fits'
        tsmapsplit['bins_per_job'] = 5
        tsmapsplit['compute_null'] = False
        tsmapsplit['outfile']      = 'cstsmapsplit_py1.dat'
        tsmapsplit['nxpix']        = 5
        tsmapsplit['nypix']        = 5
        tsmapsplit['binsz']        = 0.05
        tsmapsplit['coordsys']     = 'CEL'
        tsmapsplit['xref']         = 83.6331
        tsmapsplit['yref']         = 22.01
        tsmapsplit['proj']         = 'CAR'
        tsmapsplit['logfile']      = 'cstsmapsplit_py1.log'
        tsmapsplit['chatter']      = 2

        # Run cstsmapsplit script and save the commands
        tsmapsplit.logFileOpen()   # Make sure we get a log file
        tsmapsplit.run()
        tsmapsplit.save()

        # Check command file
        self._check_cmdfile('cstsmapsplit_py1.dat')

        # Test with higher chatter level 3
        tsmapsplit['outfile'] = 'cstsmapsplit_py2.dat'
        tsmapsplit['logfile'] = 'cstsmapsplit_py2.log'
        tsmapsplit['chatter'] = 3
        tsmapsplit.execute()

        # Check command file
        self._check_cmdfile('cstsmapsplit_py2.dat')

        # Test with higher chatter level 4
        tsmapsplit['outfile'] = 'cstsmapsplit_py3.dat'
        tsmapsplit['logfile'] = 'cstsmapsplit_py3.log'
        tsmapsplit['chatter'] = 4
        tsmapsplit.execute()

        # Check command file
        self._check_cmdfile('cstsmapsplit_py3.dat')

        # Return
        return

    # Check observation definition XML file
    def _check_cmdfile(self, filename):
        """
        Check ASCII files containing commands
        """
        # Create file name instance
        fname = gammalib.GFilename(filename)
        
        # Check if execution was successful
        self.test_assert(fname.exists(), 'Check of output filename exists')
        
        # Return
        return
