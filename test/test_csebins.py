#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cslightcrv script.
#
# Copyright (C) 2016-2017 Juergen Knoedlseder
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


# ============================= #
# Test class for csebins script #
# ============================= #
class Test(test):
    """
    Test class for csebins script

    This test class makes unit tests for the csebins script by using it from
    the command line and from Python.
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
        Set all test functions.
        """
        # Set test name
        self.name('csebins')

        # Append tests
        self.append(self._test_cmd, 'Test csebins on command line')
        self.append(self._test_python, 'Test csebins from Python')

        # Return
        return

    # Test csebins on command line
    def _test_cmd(self):
        """
        Test csebins on the command line.
        """
        # Set script name
        cslightcrv = self._script('csebins')

        # Setup cslightcrv command
        cmd = cslightcrv+' inobs="NONE"'+ \
                         ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                         ' emin=0.1 emax=100.0'+ \
                         ' outfile="csebins_cmd1.fits"'+ \
                         ' logfile="csebins_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check energy boundaries
        self._check_ebounds('csebins_cmd1.fits', 33)

        # Setup cslightcrv command
        cmd = cslightcrv+' inobs="NONE"'+ \
                         ' caldb="Database_that_does_not_exist" irf="'+self._irf+'"'+ \
                         ' emin=0.1 emax=100.0'+ \
                         ' outfile="csebins_cmd2.fits"'+ \
                         ' logfile="csebins_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csebins from Python
    def _test_python(self):
        """
        Test csebins from Python
        """
        # Set-up csebins
        ebins = cscripts.csebins()
        ebins['inobs']   = 'NONE'
        ebins['caldb']   = self._caldb
        ebins['irf']     = self._irf
        ebins['emin']    = 0.1
        ebins['emax']    = 100.0
        ebins['outfile'] = 'csebins_py1.fits'
        ebins['logfile'] = 'csebins_py1.log'
        ebins['chatter'] = 2

        # Run csebins script and energy boundaries
        ebins.logFileOpen()   # Make sure we get a log file
        ebins.run()
        ebins.save()

        # Check energy boundaries
        self._check_ebounds('csebins_py1.fits', 33)

        # Return
        return

    # Check energy boundaries result file
    def _check_ebounds(self, filename, bins):
        """
        Check energy boundaries file
        """
        # Load energy boundaries file
        ebounds = gammalib.GEbounds(filename)

        # Check number of energy boundaries
        self.test_value(ebounds.size(), bins)
        
        # Return
        return
