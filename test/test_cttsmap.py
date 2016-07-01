#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cttsmap tool.
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


# =========================== #
# Test class for cttsmap tool #
# =========================== #
class Test(test):
    """
    Test class for cttsmap tool
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
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
        self.name('cttsmap')

        # Append tests
        self.append(self._test_cmd, 'Test cttsmap on command line')
        self.append(self._test_python, 'Test cttsmap from Python')

        # Return
        return

    # Test cttsmap on command line
    def _test_cmd(self):
        """
        Test cttsmap on the command line
        """
        # Set tool name
        cttsmap = self._tool('cttsmap')

        # Setup cttsmap command
        cmd = cttsmap+' inobs="data/crab_events.fits"'+ \
                      ' inmodel="data/crab.xml" srcname="Crab"'+ \
                      ' caldb="prod2" irf="South_0.5h"'+ \
                      ' outmap="cttsmap_cmd1.fits"'+ \
                      ' nxpix=5 nypix=5 binsz=0.02'+ \
                      ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                      ' logfile="cttsmap_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('cttsmap_cmd1.fits')

        # Setup cttsmap command
        cmd = cttsmap+' inobs="event_file_that_does_not_exist.fits"'+ \
                      ' inmodel="data/crab.xml" srcname="Crab"'+ \
                      ' caldb="prod2" irf="South_0.5h"'+ \
                      ' outmap="cttsmap_cmd2.fits"'+ \
                      ' nxpix=5 nypix=5 binsz=0.02'+ \
                      ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                      ' logfile="cttsmap_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test cttsmap from Python
    def _test_python(self):
        """
        Test cttsmap from Python
        """
        # Set-up cttsmap
        tsmap = ctools.cttsmap()
        tsmap['inobs']    = 'data/crab_events.fits'
        tsmap['inmodel']  = 'data/crab.xml'
        tsmap['srcname']  = 'Crab'
        tsmap['caldb']    = 'prod2'
        tsmap['irf']      = 'South_0.5h'
        tsmap['outmap']   = 'cttsmap_py1.fits'
        tsmap['nxpix']    = 5
        tsmap['nypix']    = 5
        tsmap['binsz']    = 0.02
        tsmap['coordsys'] = 'CEL'
        tsmap['proj']     = 'CAR'
        tsmap['xref']     = 83.63
        tsmap['yref']     = 22.01
        tsmap['logfile']  = 'cttsmap_py1.log'
        tsmap['chatter']  = 2

        # Run cttsmap tool
        tsmap.logFileOpen()   # Make sure we get a log file
        tsmap.run()
        tsmap.save()

        # Check result file
        self._check_result_file('cttsmap_py1.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file
        """
        # Open result file
        fits = gammalib.GFits(filename)

        # Get HDUs
        ts        = fits['Primary']
        prefactor = fits['Prefactor']
        index     = fits['Index']

        # Check dimensions
        self.test_value(ts.naxis(), 2, 'Check for 2 dimensions')
        self.test_value(ts.naxes(0), 5, 'Check for 5 pixels in X')
        self.test_value(ts.naxes(1), 5, 'Check for 5 pixels in Y')
        self.test_value(prefactor.naxis(), 2, 'Check for 2 dimensions')
        self.test_value(prefactor.naxes(0), 5, 'Check for 5 pixels in X')
        self.test_value(prefactor.naxes(1), 5, 'Check for 5 pixels in Y')
        self.test_value(index.naxis(), 2, 'Check for 2 dimensions')
        self.test_value(index.naxes(0), 5, 'Check for 5 pixels in X')
        self.test_value(index.naxes(1), 5, 'Check for 5 pixels in Y')

        # Return
        return
