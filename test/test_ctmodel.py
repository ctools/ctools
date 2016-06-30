#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctmodel tool.
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
# Test class for ctmodel tool #
# =========================== #
class Test(test):
    """
    Test class for ctmodel tool
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
        self.name('ctmodel')

        # Append tests
        self.append(self._test_cmd, 'Test ctmodel on command line')
        self.append(self._test_python, 'Test ctmodel from Python')

        # Return
        return

    # Test ctmodel on command line
    def _test_cmd(self):
        """
        Test ctmodel on the command line
        """
        # Set tool name
        ctmodel = self._tool('ctmodel')

        # Setup ctmodel command
        cmd = ctmodel+' incube="NONE" inobs="NONE" expcube="NONE"'+\
                      ' psfcube="NONE" bkgcube="NONE"'+ \
                      ' outcube="ctmodel_cmd1.fits"'+ \
                      ' inmodel="data/crab.xml"'+ \
                      ' caldb="prod2" irf="South_0.5h"'+ \
                      ' rad=5.0 ra=83.63 dec=22.01 tmin=0.0 tmax=1800'+ \
                      ' emin=0.1 emax=100.0 enumbins=20 nxpix=200 nypix=200'+ \
                      ' binsz=0.02 coordsys="CEL" proj="CAR"'+ \
                      ' xref=83.63 yref=22.01'+ \
                      ' logfile="ctmodel_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command that does not exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctmodel_cmd1.fits')

        # Setup ctmodel command
        cmd = ctmodel+' incube="NONE" inobs="NONE" expcube="NONE"'+\
                      ' psfcube="NONE" bkgcube="NONE"'+ \
                      ' outcube="ctmodel_cmd2.fits"'+ \
                      ' inmodel="model_that_does_not_exist.xml"'+ \
                      ' caldb="prod2" irf="South_0.5h"'+ \
                      ' rad=5.0 ra=83.63 dec=22.01 tmin=0.0 tmax=1800'+ \
                      ' emin=0.1 emax=100.0 enumbins=20 nxpix=200 nypix=200'+ \
                      ' binsz=0.02 coordsys="CEL" proj="CAR"'+ \
                      ' xref=83.63 yref=22.01'+ \
                      ' logfile="ctmodel_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctmodel from Python
    def _test_python(self):
        """
        Test ctmodel from Python
        """
        # Set-up ctmodel from scratch
        model = ctools.ctmodel()
        model['incube']   = 'NONE'
        model['outcube']  = 'ctmodel_py1.fits'
        model['inmodel']  = 'data/crab.xml'
        model['inobs']    = 'NONE'
        model['expcube']  = 'NONE'
        model['psfcube']  = 'NONE'
        model['bkgcube']  = 'NONE'
        model['caldb']    = 'prod2'
        model['irf']      = 'South_0.5h'
        model['rad']      = 5
        model['ra']       = 83.63
        model['dec']      = 22.01
        model['tmin']     = 0
        model['tmax']     = 1800
        model['emin']     = 0.1
        model['emax']     = 100
        model['enumbins'] = 20
        model['nxpix']    = 200
        model['nypix']    = 200
        model['binsz']    = 0.02
        model['coordsys'] = 'CEL'
        model['proj']     = 'CAR'
        model['xref']     = 83.63
        model['yref']     = 22.01
        model['logfile']  = 'ctmodel_py1.log'
        model['chatter']  = 2

        # Run ctmodel tool
        model.logFileOpen()   # Make sure we get a log file
        model.run()
        model.save()

        # Check result file
        self._check_result_file('ctmodel_py1.fits')

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
        cube    = fits['Primary']
        ebounds = fits['EBOUNDS']
        gti     = fits['GTI']

        # Check dimensions
        self.test_value(cube.naxis(), 3, 'Check for 3 cube dimensions')
        self.test_value(cube.naxes(0), 200, 'Check for 200 pixels in X')
        self.test_value(cube.naxes(1), 200, 'Check for 200 pixels in Y')
        self.test_value(cube.naxes(2), 20, 'Check for 20 pixels in Z')

        # Return
        return
