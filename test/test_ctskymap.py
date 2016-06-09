#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctskymap tool.
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


# ============================ #
# Test class for ctskymap tool #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctskymap tool.
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
        self.name('ctskymap')

        # Append tests
        self.append(self._test_cmd, 'Test ctskymap on command line')
        self.append(self._test_python, 'Test ctskymap from Python')

        # Return
        return

    # Test ctskymap on command line
    def _test_cmd(self):
        """
        Test ctskymap on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README'):
            ctskymap = '../src/ctskymap/ctskymap'
        else:
            ctskymap = 'ctskymap'

        # Setup ctskymap command
        cmd = ctskymap+' inobs="data/crab_events.fits"'+\
                       ' outmap="ctskymap_cmd1.fits"'+ \
                       ' emin=0.1 emax=100.0 nxpix=200 nypix=200'+ \
                       ' binsz=0.02 coordsys="CEL" proj="CAR"'+ \
                       ' xref=83.63 yref=22.01'+ \
                       ' logfile="ctskymap_cmd1.log" chatter=1'

        # Execute ctskymap, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful ctskymap execution on command line')

        # Check result file
        self._check_result_file('ctskymap_cmd1.fits')

        # Setup ctskymap command
        cmd = ctskymap+' inobs="event_file_that_does_not_exist.fits"'+\
                       ' outmap="ctskymap_cmd2.fits"'+ \
                       ' emin=0.1 emax=100.0 nxpix=200 nypix=200'+ \
                       ' binsz=0.02 coordsys="CEL" proj="CAR"'+ \
                       ' xref=83.63 yref=22.01'+ \
                       ' logfile="ctskymap_cmd2.log" chatter=2'

        # Execute ctskymap, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of ctskymap execution on command line')

        # Return
        return

    # Test ctskymap from Python
    def _test_python(self):
        """
        Test ctskymap from Python.
        """
        # Set-up ctskymap
        skymap = ctools.ctskymap()
        skymap['inobs']    = 'data/crab_events.fits'
        skymap['outmap']   = 'ctskymap_py1.fits'
        skymap['emin']     = 0.1
        skymap['emax']     = 100
        skymap['nxpix']    = 200
        skymap['nypix']    = 200
        skymap['binsz']    = 0.02
        skymap['coordsys'] = 'CEL'
        skymap['proj']     = 'CAR'
        skymap['xref']     = 83.63
        skymap['yref']     = 22.01
        skymap['logfile']  = 'ctskymap_py1.log'
        skymap['chatter']  = 2

        # Run ctskymap tool
        skymap.logFileOpen()   # Make sure we get a log file
        skymap.run()
        skymap.save()

        # Check result file
        self._check_result_file('ctskymap_py1.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file.
        """
        # Open result file
        result = gammalib.GSkyMap(filename)

        # Check dimensions
        self.test_value(result.nmaps(), 1, 'Check for one map')
        self.test_value(result.nx(), 200, 'Check for 200 pixels in X')
        self.test_value(result.ny(), 200, 'Check for 200 pixels in Y')

        # Return
        return
