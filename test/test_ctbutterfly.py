#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctbutterfly tool.
#
# Copyright (C) 2014-2016 Michal Mayer
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
import ctools


# =============================== #
# Test class for ctbutterfly tool #
# =============================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctbutterfly tool.
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
        self.name('ctbutterfly')

        # Append tests
        self.append(self._test_cmd, 'Test ctbutterfly on command line')
        self.append(self._test_python, 'Test ctbutterfly from Python')

        # Return
        return

    # Test ctbutterfly on command line
    def _test_cmd(self):
        """
        Test ctbutterfly on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README.md'):
            ctbutterfly = '../src/ctbutterfly/ctbutterfly'
        else:
            ctbutterfly = 'ctbutterfly'

        # Setup ctbutterfly command
        cmd = ctbutterfly+' inobs="data/crab_events.fits"'+\
                          ' outfile="ctbutterfly_cmd1.txt"'+ \
                          ' inmodel="data/crab.xml" srcname="Crab"'+ \
                          ' caldb="prod2" irf="South_0.5h"'+ \
                          ' emin=0.1 emax=100.0'+ \
                          ' logfile="ctbutterfly_cmd1.log" chatter=1'

        # Execute ctbutterfly, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful ctbutterfly execution on command line')

        # Check result file
        self._check_result_file('ctbutterfly_cmd1.txt')

        # Setup ctbutterfly command
        cmd = ctbutterfly+' inobs="event_file_that_does_not_exist.fits"'+\
                          ' outfile="ctbutterfly_cmd2.txt"'+ \
                          ' inmodel="data/crab.xml" srcname="Crab"'+ \
                          ' caldb="prod2" irf="South_0.5h"'+ \
                          ' emin=0.1 emax=100.0'+ \
                          ' logfile="ctbutterfly_cmd2.log" chatter=1'

        # Execute ctbutterfly, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of ctbutterfly execution on command line')

        # Return
        return

    # Test ctbutterfly from Python
    def _test_python(self):
        """
        Test ctbutterfly from Python.
        """
        # Set-up ctbutterfly
        butterfly = ctools.ctbutterfly()
        butterfly['inobs']   = 'data/crab_events.fits'
        butterfly['inmodel'] = 'data/crab.xml'
        butterfly['srcname'] = 'Crab'
        butterfly['caldb']   = 'prod2'
        butterfly['irf']     = 'South_0.5h'
        butterfly['emin']    = 0.1
        butterfly['emax']    = 100.0
        butterfly['outfile'] = 'ctbutterfly_py1.txt'
        butterfly['logfile'] = 'ctbutterfly_py1.log'
        butterfly['chatter'] = 2

        # Run ctbutterfly tool
        butterfly.logFileOpen()   # Make sure we get a log file
        butterfly.run()
        butterfly.save()

        # Check result file
        self._check_result_file('ctbutterfly_py1.txt')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file.
        """
        # Open result file as CSV file
        results = gammalib.GCsv(filename)

        # Check dimensions
        self.test_value(results.nrows(), 100, 'Check for 100 rows in butterfly file')
        self.test_value(results.ncols(), 4, 'Check for 4 columns in butterfly file')

        # Return
        return
