#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comgendb script.
#
# Copyright (C) 2022 Juergen Knoedlseder
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
import comscripts
from testing import test


# ============================== #
# Test class for comgendb script #
# ============================== #
class Test(test):
    """
    Test class for comgendb script

    This test class makes unit tests for the comgendb script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test datasets and parameters
        self._archive = self._datadir + '/comptel/data'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comgendb')

        # Append tests
        self.append(self._test_cmd, 'Test comgendb on command line')
        self.append(self._test_python, 'Test comgendb from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comgendb = self._comscript('comgendb')

        # Setup command
        cmd = comgendb+' archive="'+self._archive+'" dbase="comgendb_cmd1" '+ \
                       'logfile="comgendb_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comgendb_cmd1')

        # Setup command
        cmd = comgendb+' archive="folder_that_does_not_exist" dbase="comgendb_cmd2" '+ \
                       'logfile="comgendb_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comgendb)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        db = comscripts.comgendb()
        db['archive'] = self._archive
        db['dbase']   = 'comgendb_py1'
        db['logfile'] = 'comgendb_py1.log'
        db['chatter'] = 2

        # Run script and save result
        db.logFileOpen()   # Make sure we get a log file
        db.run()
        db.save()

        # Check result
        self._check_result('comgendb_py1')

        # Return
        return

    # Check result
    def _check_result(self, dbase):
        """
        Check result
        """
        # Set database filename
        filename = '%s/dbase.fits' % (dbase)

        # Load database file
        fits = gammalib.GFits(filename)

        # Check existence of headers
        self.test_assert(fits.contains('EVP'), 'Check for "EVP" extension')
        self.test_assert(fits.contains('OAD'), 'Check for "OAD" extension')
        self.test_assert(fits.contains('TIM'), 'Check for "TIM" extension')
        self.test_assert(fits.contains('XML'), 'Check for "XML" extension')

        # Check number of rows in extension
        self.test_value(fits.table('EVP').nrows(), 1, 'Number of rows in "EVP" table')
        self.test_value(fits.table('OAD').nrows(), 1, 'Number of rows in "OAD" table')
        self.test_value(fits.table('TIM').nrows(), 0, 'Number of rows in "TIM" table')
        self.test_value(fits.table('XML').nrows(), 0, 'Number of rows in "XML" table')

        # Return
        return
