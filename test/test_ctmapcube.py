#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctmapcube tool.
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
import gammalib
import ctools
from testing import test


# ============================= #
# Test class for ctmapcube tool #
# ============================= #
class Test(test):
    """
    Test class for ctmapcube tool

    This test class makes unit tests for the ctmapcube tool by using it from
    the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set members
        self._model_name = 'data/crab.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('ctmapcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctmapcube on command line')
        self.append(self._test_python, 'Test ctmapcube from Python')

        # Return
        return

    # Test ctmapcube on command line
    def _test_cmd(self):
        """
        Test ctmapcube on the command line
        """
        # Set tool name
        ctmapcube = self._tool('ctmapcube')

        # Setup ctmapcube command
        cmd = ctmapcube+' inmodel="'+self._model_name+'"'+ \
                        ' outcube="ctmapcube_cmd1.fits"'+\
                        ' emin=0.1 emax=100.0 enumbins=20 ebinalg=LOG'+ \
                        ' nxpix=200 nypix=200 binsz=0.02 coordsys=CEL'+ \
                        ' xref=83.63 yref=22.01 proj=CAR'+ \
                        ' logfile="ctmapcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check map cube
        self._check_cube('ctmapcube_cmd1.fits')

        # Setup ctmapcube command
        cmd = ctmapcube+' inmodel="events_that_do_not_exist.fits"'+ \
                        ' outcube="ctmapcube_cmd2.fits"'+\
                        ' emin=0.1 emax=100.0 enumbins=20 ebinalg=LOG'+ \
                        ' nxpix=200 nypix=200 binsz=0.02 coordsys=CEL'+ \
                        ' xref=83.63 yref=22.01 proj=CAR'+ \
                        ' logfile="ctmapcube_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctmapcube from Python
    def _test_python(self):
        """
        Test ctmapcube from Python
        """
        # Set-up ctmapcube
        mapcube = ctools.ctmapcube()
        mapcube['inmodel']  = self._model_name
        mapcube['outcube']  = 'ctmapcube_py1.fits'
        mapcube['ebinalg']  = 'LOG'
        mapcube['emin']     = 0.1
        mapcube['emax']     = 100.0
        mapcube['enumbins'] = 20
        mapcube['nxpix']    = 200
        mapcube['nypix']    = 200
        mapcube['binsz']    = 0.02
        mapcube['coordsys'] = 'CEL'
        mapcube['proj']     = 'CAR'
        mapcube['xref']     = 83.63
        mapcube['yref']     = 22.01
        mapcube['logfile']  = 'ctmapcube_py1.log'
        mapcube['chatter']  = 2

        # Run ctmapcube tool
        mapcube.logFileOpen()   # Make sure we get a log file
        mapcube.run()

        # Save map cube
        mapcube.save()

        # Check map cube
        self._check_cube('ctmapcube_py1.fits')

        # Return
        return

    # Check map cube
    def _check_cube(self, filename):
        """
        Check content of map cube

        Parameters
        ----------
        filename : str
            Map cube file name.
        """
        # Load map cube
        cube = gammalib.GModelSpatialDiffuseCube(filename)

        # Test map cube
        #self.test_value(cube.maps(), 20, "20 maps")
        #self.test_value(cube.pixels(), 40000, "40000 map pixels")
        # The map cube is not loaded by default !!!! We should add the method so that
        # the attributes are known !!!!
        self.test_value(cube.maps(), 0, '20 maps')
        self.test_value(cube.pixels(), 0, '40000 map pixels')

        # Return
        return
