#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csadd2caldb script.
#
# Copyright (C) 2021 Juergen Knoedlseder
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
from testing import test


# ================================= #
# Test class for csadd2caldb script #
# ================================= #
class Test(test):
    """
    Test class for csadd2caldb script

    This test class makes unit tests for the csadd2caldb script by using it
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
        self._prod501 = self._datadir + '/prod5-v0.1'

        # Clear results
        os.system('rm -rf csadd2caldb_cmd1')
        os.system('rm -rf csadd2caldb_py1')
        os.system('rm -rf csadd2caldb_py1_pickle')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csadd2caldb')

        # Append tests
        self.append(self._test_cmd, 'Test csadd2caldb on command line')
        self.append(self._test_python, 'Test csadd2caldb from Python')
        self.append(self._test_pickeling, 'Test csadd2caldb pickeling')

        # Return
        return

    # Test csadd2caldb on command line
    def _test_cmd(self):
        """
        Test csadd2caldb on the command line.
        """
        # Set script name
        csadd2caldb = self._script('csadd2caldb')

        # Setup csadd2caldb command
        cmd = csadd2caldb+' indir="'+self._prod501+'"'+ \
                          ' outdir="csadd2caldb_cmd1"'+ \
                          ' logfile="csadd2caldb_cmd1.log" chatter=4'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check IRF
        self._check_irf('csadd2caldb_cmd1')

        # Setup csadd2caldb command
        cmd = csadd2caldb+' indir="folder_that_does_not_exist"'+ \
                          ' outdir="csadd2caldb_cmd2"'+ \
                          ' logfile="csadd2caldb_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check csadd2caldb --help
        self._check_help(csadd2caldb)

        # Return
        return

    # Test csadd2caldb from Python
    def _test_python(self):
        """
        Test csadd2caldb from Python
        """
        # Set-up csadd2caldb
        addirf = cscripts.csadd2caldb()
        addirf['indir']   = self._prod501
        addirf['outdir']  = 'csadd2caldb_py1'
        addirf['logfile'] = 'csadd2caldb_py1.log'

        # Run csadd2caldb script and save background model
        addirf.logFileOpen()   # Make sure we get a log file
        addirf.run()
        addirf.save()

        # Check IRF
        self._check_irf('csadd2caldb_py1')

        # Return
        return

    # Test csadd2caldb pickeling
    def _test_pickeling(self):
        """
        Test csadd2caldb pickeling
        """
        # Perform pickeling tests of empty class
        self._pickeling(cscripts.csadd2caldb())

        # Set-up csadd2caldb
        addirf = cscripts.csadd2caldb()
        addirf['indir']   = self._prod501
        addirf['outdir']  = 'csadd2caldb_py1_pickle'
        addirf['logfile'] = 'csadd2caldb_py1_pickle.log'

        # Perform pickeling tests of filled class
        obj = self._pickeling(addirf)

        # Run csadd2caldb script and save background model
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # Check IRF
        self._check_irf('csadd2caldb_py1_pickle')

        # Return
        return

    # Check IRF
    def _check_irf(self, filename, prod='prod5-v0.1', nirfs=8):
        """
        Check IRF
        """
        # Open CALDB
        caldb = gammalib.GCaldb(filename)
        caldb.open('CTA', prod)

        # Check CALDB
        self.test_value(caldb.size(), nirfs, 'Check for %d entries in CALDB' % nirfs)

        # Return
        return
