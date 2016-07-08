#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctexpcube tool.
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
# Test class for ctexpcube tool #
# ============================= #
class Test(test):
    """
    Test class for ctexpcube tool
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
        self.name('ctexpcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctexpcube on command line')
        self.append(self._test_python, 'Test ctexpcube from Python')

        # Return
        return

    # Test ctexpcube on command line
    def _test_cmd(self):
        """
        Test ctexpcube on the command line
        """
        # Set tool name
        ctexpcube = self._tool('ctexpcube')

        # Setup ctexpcube command
        cmd = ctexpcube+' inobs="data/crab_events.fits"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctexpcube_cmd1.fits"'+ \
                        ' caldb="prod2" irf="South_0.5h"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=200 nypix=200 binsz=0.02'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctexpcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctexpcube_cmd1.fits')

        # Setup ctexpcube command
        cmd = ctexpcube+' inobs="event_file_that_does_not_exist.fits"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctexpcube_cmd2.fits"'+ \
                        ' caldb="prod2" irf="South_0.5h"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=200 nypix=200 binsz=0.02'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctexpcube_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctexpcube from Python
    def _test_python(self):
        """
        Test ctexpcube from Python
        """
        # Set-up ctexpcube
        expcube = ctools.ctexpcube()
        expcube['inobs']    = 'data/crab_events.fits'
        expcube['incube']   = 'NONE'
        expcube['outcube']  = 'ctexpcube_py1.fits'
        expcube['caldb']    = 'prod2'
        expcube['irf']      = 'South_0.5h'
        expcube['ebinalg']  = 'LOG'
        expcube['emin']     = 0.1
        expcube['emax']     = 100
        expcube['enumbins'] = 20
        expcube['nxpix']    = 200
        expcube['nypix']    = 200
        expcube['binsz']    = 0.02
        expcube['coordsys'] = 'CEL'
        expcube['proj']     = 'CAR'
        expcube['xref']     = 83.63
        expcube['yref']     = 22.01
        expcube['logfile']  = 'ctexpcube_py1.log'
        expcube['chatter']  = 2

        # Run ctexpcube tool
        expcube.logFileOpen()   # Make sure we get a log file
        expcube.run()
        expcube.save()

        # Check result file
        self._check_result_file('ctexpcube_py1.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file
        """
        # Open result file
        result = gammalib.GCTACubeExposure(filename)

        # Check dimensions
        self.test_value(len(result.energies()), 21, 'Check for 21 energy values')

        # Return
        return
