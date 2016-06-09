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
import os
import gammalib
import cscripts


# ================================ #
# Test class for csfindobs script #
# ================================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for csfindobs script.

    This test class makes unit tests for the csfindobs script by using it
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
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README'):
            csfindobs = '../cscripts/csfindobs.py'
        else:
            csfindobs = 'csfindobs'

        # Setup csfindobs command
        cmd = csfindobs+' datapath="iactdata/"'+ \
                        ' prodname="unit-test"'+ \
                        ' ra=83.63 dec=22.01 rad=1.0'+ \
                        ' outfile="runlist_cmd1.dat"'+ \
                        ' logfile="csfindobs_cmd1.log" chatter=1'

        # Execute csfindobs, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful csfindobs execution on command line')

        # Check run list
        self._check_runlist('runlist_cmd1.dat')

        # Setup csfindobs command
        cmd = csfindobs+' datapath="path_that_does_not_exist"'+ \
                        ' prodname="unit-test"'+ \
                        ' ra=83.63 dec=22.01 rad=1.0'+ \
                        ' outfile="runlist_cmd2.dat"'+ \
                        ' logfile="csfindobs_cmd2.log"'

        # Execute csfindobs, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of csfindobs execution on command line')

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
