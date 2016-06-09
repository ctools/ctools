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
import os
import gammalib
import ctools


# ============================= #
# Test class for ctedispcube tool #
# ============================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctedispcube tool.
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
        self.name('ctedispcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctedispcube on command line')
        self.append(self._test_python, 'Test ctedispcube from Python')

        # Return
        return

    # Test ctedispcube on command line
    def _test_cmd(self):
        """
        Test ctedispcube on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README'):
            ctedispcube = '../src/ctedispcube/ctedispcube'
        else:
            ctedispcube = 'ctedispcube'

        # Setup ctedispcube command
        cmd = ctedispcube+' inobs="data/crab_events.fits"'+ \
                          ' incube="NONE"'+ \
                          ' outcube="ctedispcube_cmd1.fits"'+ \
                          ' caldb="prod2" irf="South_0.5h"'+ \
                          ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                          ' nxpix=10 nypix=10 binsz=0.4'+ \
                          ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                          ' migramax=2.0 migrabins=10'+ \
                          ' logfile="ctedispcube_cmd1.log" chatter=1'

        # Execute ctedispcube, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful ctedispcube execution on command line')

        # Check result file
        self._check_result_file('ctedispcube_cmd1.fits')

        # Setup ctedispcube command
        cmd = ctedispcube+' inobs="events_that_do_not_exist.fits"'+ \
                          ' incube="NONE"'+ \
                          ' outcube="ctedispcube_cmd2.fits"'+ \
                          ' caldb="prod2" irf="South_0.5h"'+ \
                          ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                          ' nxpix=10 nypix=10 binsz=0.4'+ \
                          ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                          ' migramax=2.0 migrabins=10'+ \
                          ' logfile="ctedispcube_cmd2.log" chatter=2'

        # Execute ctedispcube, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of ctedispcube execution on command line')

        # Return
        return

    # Test ctedispcube from Python
    def _test_python(self):
        """
        Test ctedispcube from Python.
        """
        # Set-up ctedispcube
        edispcube = ctools.ctedispcube()
        edispcube['inobs']     = 'data/crab_events.fits'
        edispcube['incube']    = 'NONE'
        edispcube['outcube']   = 'ctedispcube_py1.fits'
        edispcube['caldb']     = 'prod2'
        edispcube['irf']       = 'South_0.5h'
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
        Check result file.
        """
        # Open result file
        result = gammalib.GCTACubeEdisp(filename)

        # Check dimensions
        self.test_value(len(result.energies()), 21, 'Check for 21 energy maps')
        self.test_value(len(result.migras()), 10, 'Check for 10 migration bins')

        # Return
        return
