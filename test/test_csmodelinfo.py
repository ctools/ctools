#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csmodelinfo script.
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
import gammalib
import cscripts
from testing import test


# ================================= #
# Test class for csmodelinfo script #
# ================================= #
class Test(test):
    """
    Test class for csmodelinfo script

    This test class makes unit tests for the csmodelinfo script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test data
        self._models_spat = self._datadir + '/models_spatial.xml'
        self._models_spec = self._datadir + '/models_spectral.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csmodelinfo')

        # Append tests
        self.append(self._test_cmd, 'Test csmodelinfo on command line')
        self.append(self._test_python, 'Test csmodelinfo from Python')

        # Return
        return

    # Test csmodelinfo on command line
    def _test_cmd(self):
        """
        Test csmodelinfo on the command line
        """
        # Set script name
        csmodelinfo = self._script('csmodelinfo')

        # Setup csmodelinfo command
        cmd = csmodelinfo+' inmodel="'+self._models_spec+'"'+ \
                          ' outds9file="csmodelinfo_cmd1.reg"'+ \
                          ' logfile="csmodelinfo_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check region file
        self._check_ds9_file('csmodelinfo_cmd1.reg')

        # Setup csmodelinfo command
        cmd = csmodelinfo+' inmodel="model_that_does_not_exist.xml"'+ \
                          ' outds9file="csmodelinfo_cmd2.reg"'+ \
                          ' logfile="csmodelinfo_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csmodelinfo from Python
    def _test_python(self):
        """
        Test csmodelinfo from Python
        """
        # Set-up csmodelinfo
        modelinfo = cscripts.csmodelinfo()
        modelinfo['inmodel']    = self._models_spat
        modelinfo['outds9file'] = 'csmodelinfo_py1.reg'
        modelinfo['logfile']    = 'csmodelinfo_py1.log'
        modelinfo['chatter']    = 2

        # Run csmodelinfo script and save DS9 file
        modelinfo.logFileOpen()   # Make sure we get a log file
        modelinfo.run()
        modelinfo.save()

        # Check region file
        self._check_ds9_file('csmodelinfo_py1.reg')

        # Set-up csmodelinfo
        modelinfo = cscripts.csmodelinfo()
        modelinfo['inmodel']    = self._models_spec
        modelinfo['outds9file'] = 'csmodelinfo_py2.reg'
        modelinfo['logfile']    = 'csmodelinfo_py2.log'
        modelinfo['chatter']    = 3

        # Execute csmodelinfo script
        modelinfo.execute()

        # Check region file
        self._check_ds9_file('csmodelinfo_py2.reg')

        # Return
        return

    # Check region file
    def _check_ds9_file(self, filename):
        """
        Check region file
        """
        # Open file   
        f = open(filename,'r')

        # Expect "fk5" in first line
        line = f.readline().strip('\n')
        self.test_assert(line == 'fk5',
                         'Check for "fk5" in first line')

        # Expect "fk5" in first line
        line = f.readline().strip('\n')
        self.test_assert(line[0:6] == 'point(',
                         'Check for "point(" in second line')

        # Close file
        f.close()
        
        # Return
        return
