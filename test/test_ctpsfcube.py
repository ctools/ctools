#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctpsfcube tool.
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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


# ============================= #
# Test class for ctpsfcube tool #
# ============================= #
class Test(test):
    """
    Test class for ctpsfcube tool
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
        self.name('ctpsfcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctpsfcube on command line')
        self.append(self._test_python, 'Test ctpsfcube from Python')

        # Return
        return

    # Test ctpsfcube on command line
    def _test_cmd(self):
        """
        Test ctpsfcube on the command line.
        """
        # Set tool name
        ctpsfcube = self._tool('ctpsfcube')

        # Setup ctpsfcube command
        cmd = ctpsfcube+' inobs="'+self._events+'"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctpsfcube_cmd1.fits"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=200 nypix=200 binsz=0.02'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' amax=0.3 anumbins=10'+ \
                        ' logfile="ctpsfcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctpsfcube_cmd1.fits')

        # Setup ctpsfcube command
        cmd = ctpsfcube+' inobs="event_file_that_does_not_exist.fits"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctpsfcube_cmd2.fits"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=200 nypix=200 binsz=0.02'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' amax=0.3 anumbins=10'+ \
                        ' logfile="ctpsfcube_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctpsfcube from Python
    def _test_python(self):
        """
        Test ctpsfcube from Python
        """
        # Set-up ctpsfcube
        psfcube = ctools.ctpsfcube()
        psfcube['inobs']    = self._events
        psfcube['incube']   = 'NONE'
        psfcube['outcube']  = 'ctpsfcube_py1.fits'
        psfcube['caldb']    = self._caldb
        psfcube['irf']      = self._irf
        psfcube['ebinalg']  = 'LOG'
        psfcube['emin']     = 0.1
        psfcube['emax']     = 100
        psfcube['enumbins'] = 20
        psfcube['nxpix']    = 10
        psfcube['nypix']    = 10
        psfcube['binsz']    = 0.4
        psfcube['coordsys'] = 'CEL'
        psfcube['proj']     = 'CAR'
        psfcube['xref']     = 83.63
        psfcube['yref']     = 22.01
        psfcube['amax']     = 0.3
        psfcube['anumbins'] = 10
        psfcube['logfile']  = 'ctpsfcube_py1.log'
        psfcube['chatter']  = 2

        # Run ctpsfcube tool
        psfcube.logFileOpen()   # Make sure we get a log file
        psfcube.run()
        psfcube.save()

        # Check result file
        self._check_result_file('ctpsfcube_py1.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file
        """
        # Open result file
        result = gammalib.GCTACubePsf(filename)

        # Check dimensions
        self.test_value(len(result.energies()), 21, 'Check for 21 energy values')

        # Return
        return
