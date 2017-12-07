#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csmodelsois tool.
#
# Copyright (C) 2017 Josh Cardenzana
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
import cscripts
from testing import test


# =============================== #
# Test class for csmodelsois tool #
# =============================== #
class Test(test):
    """
    Test class for csmodelsois tool

    This test class makes unit tests for the csmodelsois tool by using it from
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
        self.name('csmodelsois')

        # Append tests
        self.append(self._test_cmd, 'Test csmodelsois on command line')
        self.append(self._test_python, 'Test csmodelsois from Python')

        # Return
        return

    # Test csmodelsois on command line
    def _test_cmd(self):
        """
        Test csmodelsois on the command line
        """
        # Set tool name
        csmodelsois = self._script('csmodelsois')

        # Setup csmodelsois command
        cmd = csmodelsois+' inmodel="'+self._model+'"'+ \
                          ' outcube="csmodelsois_cmd1.fits"'+\
                          ' emin=0.1 emax=100.0 enumbins=10 ebinalg=LOG'+ \
                          ' nxpix=100 nypix=100 binsz=0.04 coordsys=CEL'+ \
                          ' soilist="" outmodel=NONE'+ \
                          ' ra=83.63 dec=22.01 proj=CAR'+ \
                          ' logfile="csmodelsois_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check map cube
        self._check_result_file('csmodelsois_cmd1.fits')

        # Setup csmodelsois command
        cmd = csmodelsois+' inmodel="events_that_do_not_exist.fits"'+ \
                          ' outcube="ccsmodelsois_cmd2.fits"'+\
                          ' emin=0.1 emax=100.0 enumbins=10 ebinalg=LOG'+ \
                          ' nxpix=100 nypix=100 binsz=0.04 coordsys=CEL'+ \
                          ' soilist="" outmodel=NONE'+ \
                          ' ra=83.63 dec=22.01 proj=CAR'+ \
                          ' logfile="csmodelsois_cmd2.log" debug=yes'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Check csmodelsois --help
        self._check_help(csmodelsois)

        # Return
        return

    # Test csmodelsois from Python
    def _test_python(self):
        """
        Test csmodelsois from Python
        """
        # Allocate csmodelsois
        modelsois = cscripts.csmodelsois()

        # Check that empty csmodelsois tool holds an map cube that has no
        # energy bins
        self._check_cube(modelsois.mapcube(), nmaps=0, npixels=0)

        # Check that saving does nothing
        modelsois['outcube']  = 'csmodelsois_py0.fits'
        modelsois['outmodel'] = 'NONE'
        modelsois['logfile']  = 'csmodelsois_py0.log'
        modelsois.logFileOpen()
        modelsois.save()
        self.test_assert(not os.path.isfile('csmodelsois_py0.fits'),
             'Check that no map cube has been created')

        # Check that clearing does not lead to an exception or segfault
        modelsois.clear()

        # Now set csmodelsois parameters
        modelsois['inmodel']  = self._model
        modelsois['ebinalg']  = 'LOG'
        modelsois['emin']     = 0.1
        modelsois['emax']     = 100.0
        modelsois['enumbins'] = 10
        modelsois['nxpix']    = 100
        modelsois['nypix']    = 100
        modelsois['binsz']    = 0.04
        modelsois['coordsys'] = 'CEL'
        modelsois['proj']     = 'CAR'
        modelsois['ra']       = 83.63
        modelsois['dec']      = 22.01
        modelsois['ptsrcsig'] = 1.0
        modelsois['outcube']  = 'csmodelsois_py1.fits'
        modelsois['logfile']  = 'csmodelsois_py1.log'
        modelsois['chatter']  = 2
        modelsois['soilist']  = ''
        modelsois['outmodel'] = 'NONE'

        # Run csmodelsois tool
        modelsois.logFileOpen()   # Make sure we get a log file
        modelsois.run()

        # Save map cube
        modelsois.save()

        # Check map cube
        self._check_result_file('csmodelsois_py1.fits')

        # Allocate csmodelsois scripts, set models and do now a linear
        # binning; also do not use Gaussian point sources
        modelsois = cscripts.csmodelsois()
        modelsois['inmodel']  = self._model
        modelsois['ebinalg']  = 'LIN'
        modelsois['emin']     = 0.1
        modelsois['emax']     = 100.0
        modelsois['enumbins'] = 10
        modelsois['nxpix']    = 100
        modelsois['nypix']    = 100
        modelsois['binsz']    = 0.04
        modelsois['coordsys'] = 'CEL'
        modelsois['proj']     = 'CAR'
        modelsois['ra']       = 83.63
        modelsois['dec']      = 22.01
        modelsois['ptsrcsig'] = 0.0
        modelsois['outcube']  = 'csmodelsois_py2.fits'
        modelsois['logfile']  = 'csmodelsois_py2.log'
        modelsois['chatter']  = 4
        modelsois['soilist']  = ''
        modelsois['outmodel'] = 'NONE' 

        # Execute mapcube
        modelsois.logFileOpen()  # Needed to get a new log file
        modelsois.execute()

        # Check result file
        self._check_result_file('csmodelsois_py2.fits')

        # Update output filenames and output a new model file
        modelsois['outcube']  = 'csmodelsois_py3.fits'
        modelsois['logfile']  = 'csmodelsois_py3.log'
        modelsois['outmodel'] = 'csmodelsois_py3.xml'

        # Execute the file
        modelsois.logFileOpen()
        modelsois.execute()

        # Now check that the output model file doesnt contain a Crab Sources
        models3 = gammalib.GModels('csmodelsois_py3.xml')
        self.test_assert(not models3.contains('Crab'),
                         'Check Crab model is not present')
        self.test_assert(models3.contains(modelsois.cubemodelname()), 
                         'Check cube model is present')
        self._check_result_file('csmodelsois_py3.fits')

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
        self.test_value(len(cube.energies()), nmaps, 'Check number of energies')
        
        # Check dimensions
        self.test_value(cube.maps(), nmaps, 'Check number of maps')
        self.test_value(cube.pixels(), npixels, 'Check number of pixels')

        # Return
        return
