#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctmapcube tool.
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
import os
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
        cmd = ctmapcube+' inmodel="'+self._model+'"'+ \
                        ' outcube="ctmapcube_cmd1.fits"'+\
                        ' emin=0.1 emax=100.0 enumbins=10 ebinalg=LOG'+ \
                        ' nxpix=100 nypix=100 binsz=0.04 coordsys=CEL'+ \
                        ' xref=83.63 yref=22.01 proj=CAR'+ \
                        ' logfile="ctmapcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check map cube
        self._check_result_file('ctmapcube_cmd1.fits')

        # Setup ctmapcube command
        cmd = ctmapcube+' inmodel="events_that_do_not_exist.fits"'+ \
                        ' outcube="ctmapcube_cmd2.fits"'+\
                        ' emin=0.1 emax=100.0 enumbins=10 ebinalg=LOG'+ \
                        ' nxpix=100 nypix=100 binsz=0.04 coordsys=CEL'+ \
                        ' xref=83.63 yref=22.01 proj=CAR'+ \
                        ' logfile="ctmapcube_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctmapcube --help option
        cmd = ctmapcube+' --help'

        # Check if execution was successful in case that the CTOOLS
        # environment variable was set or failed otherwise
        if 'CTOOLS' in os.environ:
            self.test_value(self._execute(cmd), 0,
                 'Check successful execution with --help option')
        else:
            self.test_assert(self._execute(cmd) != 0,
                 'Check execution failure with --help option')

        # Return
        return

    # Test ctmapcube from Python
    def _test_python(self):
        """
        Test ctmapcube from Python
        """
        # Allocate ctmapcube
        mapcube = ctools.ctmapcube()

        # Check that empty ctmapcube tool holds an map cube that has no
        # energy bins
        self._check_cube(mapcube.mapcube(), nmaps=0, npixels=0)

        # Check that saving does not nothing
        mapcube['outcube'] = 'ctmapcube_py0.fits'
        mapcube['logfile'] = 'ctmapcube_py0.log'
        mapcube.logFileOpen()
        mapcube.save()
        self.test_assert(not os.path.isfile('ctmapcube_py0.fits'),
             'Check that no map cube has been created')

        # Check that clearing does not lead to an exception or segfault
        mapcube.clear()

        # Now set ctmapcube parameters
        mapcube['inmodel']  = self._model
        mapcube['ebinalg']  = 'LOG'
        mapcube['emin']     = 0.1
        mapcube['emax']     = 100.0
        mapcube['enumbins'] = 10
        mapcube['nxpix']    = 100
        mapcube['nypix']    = 100
        mapcube['binsz']    = 0.04
        mapcube['coordsys'] = 'CEL'
        mapcube['proj']     = 'CAR'
        mapcube['xref']     = 83.63
        mapcube['yref']     = 22.01
        mapcube['outcube']  = 'ctmapcube_py1.fits'
        mapcube['logfile']  = 'ctmapcube_py1.log'
        mapcube['chatter']  = 2

        # Run ctmapcube tool
        mapcube.logFileOpen()   # Make sure we get a log file
        mapcube.run()

        # Save map cube
        mapcube.save()

        # Check map cube
        self._check_result_file('ctmapcube_py1.fits')

        # Copy ctmapcube tool
        cpy_mapcube = mapcube.copy()

        # Check map cube of ctmapcube copy
        self._check_cube(cpy_mapcube.mapcube())

        # Execute copy of ctmapcube tool again, now with a higher chatter
        # level than before
        cpy_mapcube['outcube'] = 'ctmapcube_py2.fits'
        cpy_mapcube['logfile'] = 'ctmapcube_py2.log'
        cpy_mapcube['publish'] = True
        cpy_mapcube['chatter'] = 3
        cpy_mapcube.logFileOpen()  # Needed to get a new log file
        cpy_mapcube.execute()

        # Check result file
        self._check_result_file('ctmapcube_py2.fits')

        # Now clear copy of ctmapcube tool
        cpy_mapcube.clear()

        # Check that the cleared copy has also cleared the map cube
        self._check_cube(cpy_mapcube.mapcube(), nmaps=0, npixels=0)

        # Allocate ctmapcube tools, set models and do now a linear
        # binning; also do not use Gaussian point sources
        mapcube = ctools.ctmapcube()
        mapcube.models(gammalib.GModels(self._model))
        mapcube['ebinalg']  = 'LIN'
        mapcube['emin']     = 0.1
        mapcube['emax']     = 100.0
        mapcube['enumbins'] = 10
        mapcube['nxpix']    = 100
        mapcube['nypix']    = 100
        mapcube['binsz']    = 0.04
        mapcube['coordsys'] = 'CEL'
        mapcube['proj']     = 'CAR'
        mapcube['xref']     = 83.63
        mapcube['yref']     = 22.01
        mapcube['ptsrcsig'] = 0.0
        mapcube['outcube']  = 'ctmapcube_py3.fits'
        mapcube['logfile']  = 'ctmapcube_py3.log'
        mapcube['chatter']  = 4

        # Execute mapcube
        mapcube.logFileOpen()  # Needed to get a new log file
        mapcube.execute()

        # Check result file
        self._check_result_file('ctmapcube_py3.fits')

        # Publish with a name
        mapcube.publish('My map cube')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check content of map cube

        Parameters
        ----------
        filename : str
            Map cube file name
        """
        # Load map cube
        cube = gammalib.GModelSpatialDiffuseCube(filename)

        # Check map cube
        self._check_cube(cube)

        # Return
        return

    # Check map cube
    def _check_cube(self, cube, nmaps=11, npixels=10000):
        """
        Check map cube

        Parameters
        ----------
        cube : `~gammalib.GModelSpatialDiffuseCube`
            Map cube
        nmaps : int, optional
            Number of maps
        npixels : int, optional
            Number of pixels
        """
        # Get energies (this forces loading in case the map cube is not
        # loaded)
        energies = cube.energies()
        self.test_value(len(cube.energies()), nmaps, 'Check number of energies')
        
        # Check dimensions
        self.test_value(cube.maps(), nmaps, 'Check number of maps')
        self.test_value(cube.pixels(), npixels, 'Check number of pixels')

        # Return
        return
