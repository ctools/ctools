#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctbkgcube tool.
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
import os
import gammalib
import ctools


# ============================= #
# Test class for ctbkgcube tool #
# ============================= #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctbkgcube tool.
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
        self.name('ctbkgcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctbkgcube on command line')
        self.append(self._test_python, 'Test ctbkgcube from Python')

        # Return
        return

    # Test ctbkgcube on command line
    def _test_cmd(self):
        """
        Test ctbkgcube on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README'):
            ctbkgcube = '../src/ctbkgcube/ctbkgcube'
        else:
            ctbkgcube = 'ctbkgcube'

        # Setup ctbkgcube command
        cmd = ctbkgcube+' inobs="data/crab_events.fits"'+ \
                        ' inmodel="data/crab.xml"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctbkgcube_cmd1.fits"'+ \
                        ' outmodel="ctbkgcube_cmd1.xml"'+ \
                        ' caldb="prod2" irf="South_0.5h"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=10 nypix=10 binsz=0.4'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctbkgcube_cmd1.log" chatter=1'

        # Execute ctbkgcube, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful ctbkgcube execution on command line')

        # Check result file
        self._check_result_file('ctbkgcube_cmd1.fits')

        # Setup ctbkgcube command
        cmd = ctbkgcube+' inobs="event_file_that_does_not_exist.fits"'+ \
                        ' inmodel="data/crab.xml"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctbkgcube_cmd2.fits"'+ \
                        ' outmodel="ctbkgcube_cmd2.xml"'+ \
                        ' caldb="prod2" irf="South_0.5h"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=10 nypix=10 binsz=0.4'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctbkgcube_cmd2.log" chatter=1'

        # Execute ctbkgcube, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of ctbkgcube execution on command line')

        # Return
        return

    # Test ctbkgcube from Python
    def _test_python(self):
        """
        Test ctbkgcube from Python.
        """
        # Set-up ctbkgcube
        bkgcube = ctools.ctbkgcube()
        bkgcube['inobs']    = 'data/crab_events.fits'
        bkgcube['inmodel']  = 'data/crab.xml'
        bkgcube['incube']   = 'NONE'
        bkgcube['outcube']  = 'ctbkgcube_py1.fits'
        bkgcube['outmodel'] = 'ctbkgcube_py1.xml'
        bkgcube['caldb']    = 'prod2'
        bkgcube['irf']      = 'South_0.5h'
        bkgcube['ebinalg']  = 'LOG'
        bkgcube['emin']     = 0.1
        bkgcube['emax']     = 100
        bkgcube['enumbins'] = 20
        bkgcube['nxpix']    = 10
        bkgcube['nypix']    = 10
        bkgcube['binsz']    = 0.4
        bkgcube['coordsys'] = 'CEL'
        bkgcube['proj']     = 'CAR'
        bkgcube['xref']     = 83.63
        bkgcube['yref']     = 22.01
        bkgcube['logfile']  = 'ctbkgcube_py1.log'
        bkgcube['chatter']  = 2

        # Run ctbkgcube tool
        bkgcube.logFileOpen()   # Make sure we get a log file
        bkgcube.run()
        bkgcube.save()

        # Check result file
        self._check_result_file('ctbkgcube_py1.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file.
        """
        # Open result file
        result = gammalib.GCTACubeBackground(filename)

        # Check dimensions
        self.test_value(len(result.energies()), 21, 'Check for 21 energy values')

        # Return
        return
