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
        self.name('csfindobs')

        # Append tests
        self.append(self._test_cmd, 'Test csfindobs on command line')
        self.append(self._test_python, 'Test csfindobs from Python')

        # Return
        return

    # Test csfindobs on command line
    def _test_cmd(self):
        """
        Test csfindobs on the command line
        """
        # Set script name
        csfindobs = self._script('csfindobs')

        # Setup csfindobs command
        cmd = csfindobs+' datapath="'+self._datapath+'"'+ \
                        ' prodname="unit-test"'+ \
                        ' ra=83.63 dec=22.01 rad=1.5'+ \
                        ' outfile="runlist_cmd1.dat"'+ \
                        ' logfile="csfindobs_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check run list
        self._check_runlist('runlist_cmd1.dat',5)

        # Setup csfindobs command
        cmd = csfindobs+' datapath="path_that_does_not_exist"'+ \
                        ' prodname="unit-test"'+ \
                        ' ra=83.63 dec=22.01 rad=1.5'+ \
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
        Test csfindobs from Python
        """
        # Set-up csfindobs
        findobs = cscripts.csfindobs()
        findobs['datapath'] = self._datapath
        findobs['prodname'] = 'unit-test'
        findobs['ra']       = 83.63
        findobs['dec']      = 22.01
        findobs['rad']      = 1.5
        findobs['outfile']  = 'runlist_py1.dat'
        findobs['logfile']  = 'csfindobs_py1.log'
        findobs['chatter']  = 4

        # Run csfindobs script
        findobs.logFileOpen()   # Make sure we get a log file
        findobs.run()

        # Check list of runs
        runs = findobs.runs()
        self.test_value(len(runs), 5,
                        'Check for number of runs in runlist')
        self.test_value(str(runs[0]), '15000',
                        'Check run ID in runlist')

        # Save run list
        findobs.save()

        # Check run list
        self._check_runlist('runlist_py1.dat', runs=5)

        # Set-up csfindobs for undefined spatial parameters
        findobs = cscripts.csfindobs()
        findobs['datapath'] = self._datapath
        findobs['prodname'] = 'unit-test'
        findobs['ra']       = 'UNDEFINED'
        findobs['dec']      = 'UNDEFINED'
        findobs['rad']      = 'UNDEFINED'
        findobs['outfile']  = 'runlist_py2.dat'
        findobs['logfile']  = 'csfindobs_py2.log'
        findobs['chatter']  = 3

        # Run csfindobs script and save run list
        findobs.execute()

        # Check run list
        self._check_runlist('runlist_py2.dat', runs=9)

        # Set-up csfindobs for user expression
        findobs = cscripts.csfindobs()
        findobs['datapath']   = self._datapath
        findobs['prodname']   = 'unit-test'
        findobs['ra']         = 329.71
        findobs['dec']        = -30.2339
        findobs['rad']        = 2.0
        findobs['expression'] = 'ONTIME<5.0'
        findobs['outfile']    = 'runlist_py3.dat'
        findobs['logfile']    = 'csfindobs_py3.log'
        findobs['chatter']    = 4

        # Run csfindobs script and save run list
        findobs.execute()

        # Check run list
        self._check_runlist('runlist_py3.dat', runs=0)

        # Return
        return

    # Check run list result file
    def _check_runlist(self, filename, runs=1):
        """
        Check run list file

        Parameters
        ----------
        filename : str
            Run list file name
        runs : int, optional
            Expected number of runs in runlist file
        """
        # Open run list file as CSV file
        runlist = gammalib.GCsv(filename)

        # Check dimensions
        self.test_value(runlist.nrows(), runs,
                        'Check for number of runs in runlist file')
        if runs > 0:
            self.test_value(runlist.ncols(), 1,
                            'Check for single column in runlist file')
            self.test_value(runlist[0,0], '15000',
                            'Check for run ID in runlist file')
        
        # Return
        return
