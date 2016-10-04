#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csroot2caldb script.
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
from testing import test


# ================================== #
# Test class for csroot2caldb script #
# ================================== #
class Test(test):
    """
    Test class for csroot2caldb script

    This test class makes unit tests for the cslightcrv script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set members
        self._infile = self._datadir + '/irf.root'

        # Clean calibration database
        os.system('rm -rf csroot2caldb_* >/dev/null 2>&1')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csroot2caldb')

        # Append tests
        self.append(self._test_cmd, 'Test csroot2caldb on command line')
        self.append(self._test_python, 'Test csroot2caldb from Python')

        # Return
        return

    # Test csroot2caldb on command line
    def _test_cmd(self):
        """
        Test csroot2caldb on the command line
        """
        # Set script name
        csroot2caldb = self._script('csroot2caldb')

        # Setup csroot2caldb command
        cmd = csroot2caldb+' infile="'+self._infile+'"'+ \
                           ' outdir="csroot2caldb_cmd1"'+ \
                           ' inst=prod3 id=South_50h version=prod3'+ \
                           ' analysis=DESY zenith=20.0 azimuth=0.0'+ \
                           ' logfile="csroot2caldb_cmd1.log" chatter=2'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Return
        return

    # Test csroot2caldb from Python
    def _test_python(self):
        """
        Test csroot2caldb from Python
        """
        # Set-up csroot2caldb with a maximum of parameters so that the
        # corresponding code gets tests
        caldb = cscripts.csroot2caldb()
        caldb['infile']        = self._infile
        caldb['outdir']        = 'csroot2caldb_py1'
        caldb['inst']          = 'prod3'
        caldb['id']            = 'South_50h'
        caldb['version']       = 'prod3'
        caldb['analysis']      = 'DESY'
        caldb['zenith']        = 20.0
        caldb['azimuth']       = 0.0
        caldb['psftype']       = 'King'
        caldb['norm1d']        = True
        caldb['rebin']         = True
        caldb['eascale']       = 1.1
        caldb['bgdscale']      = 0.9
        caldb['bgdoversample'] = 3
        caldb['bgdethres']     = 10.0
        caldb['bgdinfill']     = True
        caldb['logfile']       = 'csroot2caldb_py1.log'
        caldb['chatter']       = 4

        # Run script
        caldb.logFileOpen()   # Make sure we get a log file
        caldb.run()

        # Return
        return
