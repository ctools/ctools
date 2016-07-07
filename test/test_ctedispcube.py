#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctedispcube tool.
#
# Copyright (C) 2016 Maria Haupt
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
import ctools
from testing import test


# =============================== #
# Test class for ctedispcube tool #
# =============================== #
class Test(test):
    """
    Test class for ctedispcube tool
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
        self.name('ctedispcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctedispcube on command line')
        self.append(self._test_python, 'Test ctedispcube from Python')

        # Return
        return

    # Test ctedispcube on command line
    def _test_cmd(self):
        """
        Test ctedispcube on the command line
        """
        # Set tool name
        ctedispcube = self._tool('ctedispcube')

        # Setup ctedispcube command
        cmd = ctedispcube+' inobs="'+self._events+'"'+ \
                          ' incube="NONE"'+ \
                          ' outcube="ctedispcube_cmd1.fits"'+ \
                          ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                          ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                          ' nxpix=10 nypix=10 binsz=0.4'+ \
                          ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                          ' migramax=2.0 migrabins=10'+ \
                          ' logfile="ctedispcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctedispcube_cmd1.fits')

        # Setup ctedispcube command
        cmd = ctedispcube+' inobs="events_that_do_not_exist.fits"'+ \
                          ' incube="NONE"'+ \
                          ' outcube="ctedispcube_cmd2.fits"'+ \
                          ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                          ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                          ' nxpix=10 nypix=10 binsz=0.4'+ \
                          ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                          ' migramax=2.0 migrabins=10'+ \
                          ' logfile="ctedispcube_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctedispcube from Python
    def _test_python(self):
        """
        Test ctedispcube from Python
        """
        # Set-up ctedispcube
        edispcube = ctools.ctedispcube()
        edispcube['inobs']     = self._events
        edispcube['incube']    = 'NONE'
        edispcube['outcube']   = 'ctedispcube_py1.fits'
        edispcube['caldb']     = self._caldb
        edispcube['irf']       = self._irf
        edispcube['ebinalg']   = 'LOG'
        edispcube['emin']      = 0.1
        edispcube['emax']      = 100
        edispcube['enumbins']  = 20
        edispcube['nxpix']     = 10
        edispcube['nypix']     = 10
        edispcube['binsz']     = 0.4
        edispcube['coordsys']  = 'CEL'
        edispcube['proj']      = 'CAR'
        edispcube['xref']      = 83.63
        edispcube['yref']      = 22.01
        edispcube['migramax']  = 2.0
        edispcube['migrabins'] = 10
        edispcube['logfile']   = 'ctedispcube_py1.log'
        edispcube['chatter']   = 2

        # Run ctedispcube tool
        edispcube.logFileOpen()   # Make sure we get a log file
        edispcube.run()
        edispcube.save()

        # Check result file
        self._check_result_file('ctedispcube_py1.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file
        """
        # Open result file
        result = gammalib.GCTACubeEdisp(filename)

        # Check dimensions
        self.test_value(len(result.energies()), 21, 'Check for 21 energy maps')
        self.test_value(len(result.migras()), 10, 'Check for 10 migration bins')

        # Return
        return
