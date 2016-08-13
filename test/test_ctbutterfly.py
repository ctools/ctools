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
from testing import test


# =============================== #
# Test class for ctbutterfly tool #
# =============================== #
class Test(test):
    """
    Test class for ctbutterfly tool
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
        self.name('ctbutterfly')

        # Append tests
        self.append(self._test_cmd, 'Test ctbutterfly on command line')
        self.append(self._test_python, 'Test ctbutterfly from Python')

        # Return
        return

    # Test ctbutterfly on command line
    def _test_cmd(self):
        """
        Test ctbutterfly on the command line
        """
        # Set tool name
        ctbutterfly = self._tool('ctbutterfly')

        # Setup ctbutterfly command
        cmd = ctbutterfly+' inobs="'+self._events+'"'+\
                          ' outfile="ctbutterfly_cmd1.dat"'+ \
                          ' inmodel="'+self._model +'" srcname="Crab"'+ \
                          ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                          ' emin=0.1 emax=100.0'+ \
                          ' logfile="ctbutterfly_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctbutterfly_cmd1.dat')

        # Setup ctbutterfly command
        cmd = ctbutterfly+' inobs="event_file_that_does_not_exist.fits"'+\
                          ' outfile="ctbutterfly_cmd2.dat"'+ \
                          ' inmodel="'+self._model +'" srcname="Crab"'+ \
                          ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                          ' emin=0.1 emax=100.0'+ \
                          ' logfile="ctbutterfly_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctbutterfly --help option
        cmd = ctbutterfly+' --help'

        # Check if execution was successful in case that the CTOOLS
        # environment variable was set or failed otherwise
        if 'CTOOLS' in os.environ:
            self.test_value(self._execute(cmd), 0,
                 'Check successful execution with --help option')
        else:
            self.test_assert(self._execute(cmd) != 0,
                 'Check execution failure with --help option')

        # Return
        return

    # Test ctbutterfly from Python
    def _test_python(self):
        """
        Test ctbutterfly from Python
        """
        # Set-up ctbutterfly
        butterfly = ctools.ctbutterfly()
        butterfly['inobs']   = self._events
        butterfly['inmodel'] = self._model
        butterfly['srcname'] = 'Crab'
        butterfly['caldb']   = self._caldb
        butterfly['irf']     = self._irf
        butterfly['emin']    = 0.1
        butterfly['emax']    = 100.0
        butterfly['outfile'] = 'ctbutterfly_py1.dat'
        butterfly['logfile'] = 'ctbutterfly_py1.log'
        butterfly['chatter'] = 2

        # Run ctbutterfly tool
        butterfly.logFileOpen()   # Make sure we get a log file
        butterfly.run()
        butterfly.save()

        # Check result file
        self._check_result_file('ctbutterfly_py1.dat')

        # Set-up ctbutterfly
        butterfly = ctools.ctbutterfly()
        butterfly['inobs']   = self._events
        butterfly['inmodel'] = self._model
        butterfly['srcname'] = 'Crab'
        butterfly['caldb']   = self._caldb
        butterfly['irf']     = self._irf
        butterfly['emin']    = 0.1
        butterfly['emax']    = 100.0
        butterfly['fit']     = True
        butterfly['outfile'] = 'ctbutterfly_py2.dat'
        butterfly['logfile'] = 'ctbutterfly_py2.log'
        butterfly['chatter'] = 3

        # Execute ctbutterfly tool
        butterfly.logFileOpen()   # Make sure we get a log file
        butterfly.execute()

        # Check result file
        self._check_result_file('ctbutterfly_py2.dat')

        # Recover observation container
        obs = butterfly.obs()

        # TODO: Do someting test on observation container

        # Copy ctbkgcube tool and execute copy
        cpy_butterfly = butterfly
        cpy_butterfly['outfile'] = 'ctbutterfly_py3.dat'
        cpy_butterfly['logfile'] = 'ctbutterfly_py3.log'
        cpy_butterfly['chatter'] = 4
        cpy_butterfly.execute()

        # Check result file
        self._check_result_file('ctbutterfly_py3.dat')

        # Clear ctbkgcube tool
        butterfly.clear()

        # TODO: Do some test after clearing

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file
        """
        # Open result file as CSV file
        results = gammalib.GCsv(filename)

        # Check dimensions
        self.test_value(results.nrows(), 100,
             'Check for 100 rows in butterfly file')
        self.test_value(results.ncols(), 4,
             'Check for 4 columns in butterfly file')

        # Return
        return
