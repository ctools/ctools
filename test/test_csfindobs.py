#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csfindobs script.
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


# =============================== #
# Test class for csfindobs script #
# =============================== #
class Test(test):
    """
    Test class for csfindobs script

    This test class makes unit tests for the csfindobs script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        test.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csfindobs')

        # Append tests
        self.append(self._test_cmd, 'Test csfindobs on command line')
        self.append(self._test_python, 'Test csfindobs from Python')

        # Return
        return

    # Test csfindobs on command line
    def _test_cmd(self):
        """
        Test csfindobs on the command line.
        """
        # Set script name
        csfindobs = self._script('csfindobs')

        # Setup csfindobs command
        cmd = csfindobs+' datapath="iactdata/"'+ \
                        ' prodname="unit-test"'+ \
                        ' ra=83.63 dec=22.01 rad=1.0'+ \
                        ' outfile="runlist_cmd1.dat"'+ \
                        ' logfile="csfindobs_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command that does not exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check run list
        self._check_runlist('runlist_cmd1.dat')

        # Setup csfindobs command
        cmd = csfindobs+' datapath="path_that_does_not_exist"'+ \
                        ' prodname="unit-test"'+ \
                        ' ra=83.63 dec=22.01 rad=1.0'+ \
                        ' outfile="runlist_cmd2.dat"'+ \
                        ' logfile="csfindobs_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csfindobs from Python
    def _test_python(self):
        """
        Test csfindobs from Python.
        """
        # Set-up csfindobs
        findobs = cscripts.csfindobs()
        findobs['datapath'] = 'iactdata/'
        findobs['prodname'] = 'unit-test'
        findobs['outfile']  = 'runlist_py1.dat'
        findobs['ra']       = 83.63
        findobs['dec']      = 22.01
        findobs['rad']      = 1.0
        findobs['logfile']  = 'csfindobs_py1.log'
        findobs['chatter']  = 2

        # Run csfindobs script and save run list
        findobs.logFileOpen()   # Make sure we get a log file
        findobs.run()
        findobs.save()

        # Check run list
        self._check_runlist('runlist_py1.dat')

        # Return
        return

    # Check run list result file
    def _check_runlist(self, filename):
        """
        Check run list file.
        """
        # Open run list file as CSV file
        runlist = gammalib.GCsv(filename)

        # Check dimensions
        self.test_value(runlist.ncols(), 1,
                        'Check for single column in runlist file')
        self.test_value(runlist.nrows(), 1,
                        'Check for single row in runlist file')
        self.test_assert(runlist[0,0] == '0',
                         'Check for value "0" in runlist file')
        
        # Return
        return
