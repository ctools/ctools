#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csiactdata script
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


# ================================ #
# Test class for csiactdata script #
# ================================ #
class Test(test):
    """
    Test class for csiactdata script

    This test class makes unit tests for the csiactdata script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        test.__init__(self)

        # Set data members
        self._datapath = self._datadir + '/iactdata'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csiactdata')

        # Append tests
        self.append(self._test_cmd, 'Test csiactdata on command line')
        self.append(self._test_python, 'Test csiactdata from Python')

        # Return
        return

    # Test csiactdata on command line
    def _test_cmd(self):
        """
        Test csiactdata on the command line
        """
        # Set script name
        csiactdata = self._script('csiactdata')

        # Setup csiactdata command
        cmd = csiactdata+' datapath="'+self._datapath+'"'+ \
                         ' logfile="csiactdata_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Setup csiactdata command
        cmd = csiactdata+' datapath="data_path_that_does_not_exist"'+ \
                         ' logfile="csiactdata_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csiactdata from Python
    def _test_python(self):
        """
        Test csiactdata from Python
        """
        # Set-up csiactdata
        iactdata = cscripts.csiactdata()
        iactdata['datapath'] = self._datapath
        iactdata['logfile']  = 'csiactdata_py1.log'
        iactdata['chatter']  = 2

        # Run csiactdata script and save run list
        iactdata.logFileOpen()   # Make sure we get a log file
        iactdata.run()
        
        # Get available configs
        names = iactdata.names()
        
        # Check if there is at least one config available
        self.test_assert(len(names) > 0,
             'Check successful execution from Python')

        # Set-up csiactdata
        iactdata = cscripts.csiactdata()
        iactdata['datapath'] = self._datapath
        iactdata['logfile']  = 'csiactdata_py2.log'
        iactdata['chatter']  = 4

        # Run csiactdata script and save run list
        iactdata.execute()
        
        # Get available configs
        names = iactdata.names()
        
        # Check if there is at least one config available
        self.test_assert(len(names) > 0,
             'Check successful execution from Python')

        # Return
        return
