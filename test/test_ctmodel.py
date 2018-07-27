#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctmodel tool.
#
# Copyright (C) 2014-2018 Juergen Knoedlseder
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

        # Set test data
        self._obs_unbinned = self._datadir + '/obs_unbinned.xml'
        self._expcube      = self._datadir + '/crab_expcube.fits'
        self._psfcube      = self._datadir + '/crab_psfcube.fits'
        self._edispcube    = self._datadir + '/crab_edispcube.fits'
        self._bkgcube      = self._datadir + '/crab_bkgcube.fits'
        self._bkgmodel     = self._datadir + '/crab_bkgcube.xml'

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
                      ' inmodel="'+self._model+'"'+ \
                      ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                      ' rad=5.0 ra=83.63 dec=22.01 tmin=0.0 tmax=1800'+ \
                      ' emin=1.0 emax=100.0 enumbins=5  ebinalg="LOG"'+ \
                      ' nxpix=40 nypix=40'+ \
                      ' binsz=0.1 coordsys="CEL" proj="CAR"'+ \
                      ' xref=83.63 yref=22.01'+ \
                      ' logfile="ctmodel_cmd1.log" chatter=1'

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
                      ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                      ' rad=5.0 ra=83.63 dec=22.01 tmin=0.0 tmax=1800'+ \
                      ' emin=1.0 emax=100.0 enumbins=5  ebinalg="LOG"'+ \
                      ' nxpix=40 nypix=40'+ \
                      ' binsz=0.1 coordsys="CEL" proj="CAR"'+ \
                      ' xref=83.63 yref=22.01'+ \
                      ' logfile="ctmodel_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check ctmodel --help
        self._check_help(ctmodel)

        # Return
        return

    # Test ctmodel from Python
    def _test_python(self):
        """
        Test ctmodel from Python
        """
        # Allocate ctmodel
        model = ctools.ctmodel()

        # Check that empty ctmodel tool holds an empty model cube
        self._check_cube(model.cube(), nx=0, ny=0, nebins=0)

        # Check that saving does not nothing
        model['outcube']  = 'ctmodel_py0.fits'
        model['logfile']  = 'ctmodel_py0.log'
        model.logFileOpen()
        model.save()
        self.test_assert(not os.path.isfile('ctmodel_py0.fits'),
             'Check that no model cube has been created')

        # Check that clearing does not lead to an exception or segfault
        model.clear()

        # Now set ctmodel parameters
        model['incube']   = 'NONE'
        model['inmodel']  = self._model
        model['inobs']    = 'NONE'
        model['expcube']  = 'NONE'
        model['psfcube']  = 'NONE'
        model['bkgcube']  = 'NONE'
        model['caldb']    = self._caldb
        model['irf']      = self._irf
        model['rad']      = 5
        model['ra']       = 83.63
        model['dec']      = 22.01
        model['tmin']     = 0
        model['tmax']     = 1800
        model['ebinalg']  = 'LOG'
        model['emin']     = 1.0
        model['emax']     = 100.0
        model['enumbins'] = 5
        model['nxpix']    = 40
        model['nypix']    = 40
        model['binsz']    = 0.1
        model['coordsys'] = 'CEL'
        model['proj']     = 'CAR'
        model['xref']     = 83.63
        model['yref']     = 22.01
        model['outcube']  = 'ctmodel_py1.fits'
        model['logfile']  = 'ctmodel_py1.log'
        model['chatter']  = 2

        # Run ctmodel tool
        model.logFileOpen()   # Make sure we get a log file
        model.run()
        model.save()

        # Check result file
        self._check_result_file('ctmodel_py1.fits')

        # Copy ctmodel tool
        cpy_model = model.copy()

        # Check model cube of ctmodel copy
        self._check_cube(cpy_model.cube())

        # Execute copy of ctmodel tool again, now with a higher chatter
        # level than before
        cpy_model['outcube'] = 'ctmodel_py2.fits'
        cpy_model['logfile'] = 'ctmodel_py2.log'
        cpy_model['publish'] = True
        cpy_model['chatter'] = 3
        cpy_model.logFileOpen()  # Needed to get a new log file
        cpy_model.execute()

        # Check result file
        self._check_result_file('ctmodel_py2.fits')

        # Now clear copy of ctmodel tool
        cpy_model.clear()

        # Check that the cleared copy has also cleared the model cube
        self._check_cube(cpy_model.cube(), nx=0, ny=0, nebins=0)

        # Get mixed observation container
        obs = self._obs_mixed()

        # Setup ctmodel tool from observation container and get model cube
        # definition from counts cube
        model = ctools.ctmodel(obs)
        model.models(gammalib.GModels(self._model))
        model['incube']   = self._cntcube
        model['expcube']  = 'NONE'
        model['psfcube']  = 'NONE'
        model['bkgcube']  = 'NONE'
        model['caldb']    = self._caldb
        model['irf']      = self._irf
        model['outcube']  = 'ctmodel_py3.fits'
        model['logfile']  = 'ctmodel_py3.log'
        model['chatter']  = 4

        # Execute ctmodel tool
        model.logFileOpen()  # Needed to get a new log file
        model.execute()

        # Check result file
        self._check_result_file('ctmodel_py3.fits', nx=20, ny=20, nebins=5)

        # Publish with name
        model.publish('My model')

        # Setup ctmodel from event list FITS file
        model = ctools.ctmodel()
        model['inobs']     = self._events
        model['incube']    = 'NONE'
        model['inmodel']   = self._model
        model['expcube']   = 'NONE'
        model['psfcube']   = 'NONE'
        model['edispcube'] = 'NONE'
        model['bkgcube']   = 'NONE'
        model['caldb']     = self._caldb
        model['irf']       = self._irf
        model['rad']       = 5
        model['ra']        = 83.63
        model['dec']       = 22.01
        model['tmin']      = 0
        model['tmax']      = 1800
        model['ebinalg']   = 'LOG'
        model['emin']      = 1.0
        model['emax']      = 100.0
        model['enumbins']  = 5
        model['nxpix']     = 40
        model['nypix']     = 40
        model['binsz']     = 0.1
        model['coordsys']  = 'CEL'
        model['proj']      = 'CAR'
        model['xref']      = 83.63
        model['yref']      = 22.01
        model['outcube']   = 'ctmodel_py4.fits'
        model['logfile']   = 'ctmodel_py4.log'
        model['chatter']   = 2

        # Run ctmodel tool
        model.logFileOpen()   # Make sure we get a log file
        model.execute()

        # Check result file
        self._check_result_file('ctmodel_py4.fits')

        # Setup ctmodel from observation definition XML file
        model = ctools.ctmodel()
        model['inobs']     = self._obs_unbinned
        model['incube']    = 'NONE'
        model['inmodel']   = self._model
        model['expcube']   = 'NONE'
        model['psfcube']   = 'NONE'
        model['edispcube'] = 'NONE'
        model['bkgcube']   = 'NONE'
        model['caldb']     = self._caldb
        model['irf']       = self._irf
        model['rad']       = 5
        model['ra']        = 83.63
        model['dec']       = 22.01
        model['tmin']      = 0
        model['tmax']      = 1800
        model['ebinalg']   = 'LOG'
        model['emin']      = 1.0
        model['emax']      = 100.0
        model['enumbins']  = 5
        model['nxpix']     = 40
        model['nypix']     = 40
        model['binsz']     = 0.1
        model['coordsys']  = 'CEL'
        model['proj']      = 'CAR'
        model['xref']      = 83.63
        model['yref']      = 22.01
        model['outcube']   = 'ctmodel_py5.fits'
        model['logfile']   = 'ctmodel_py5.log'
        model['chatter']   = 2

        # Run ctmodel tool
        model.logFileOpen()   # Make sure we get a log file
        model.execute()

        # Check result file
        self._check_result_file('ctmodel_py5.fits')

        # Setup ctmodel for stacked IRF
        model = ctools.ctmodel()
        model['inobs']     = 'NONE'
        model['incube']    = 'NONE'
        model['inmodel']   = self._bkgmodel
        model['expcube']   = self._expcube
        model['psfcube']   = self._psfcube
        model['edispcube'] = 'NONE'
        model['bkgcube']   = self._bkgcube
        model['edisp']     = False
        model['rad']       = 5
        model['ra']        = 83.63
        model['dec']       = 22.01
        model['tmin']      = 0
        model['tmax']      = 1800
        model['ebinalg']   = 'LOG'
        model['emin']      = 1.0
        model['emax']      = 100.0
        model['enumbins']  = 5
        model['nxpix']     = 40
        model['nypix']     = 40
        model['binsz']     = 0.1
        model['coordsys']  = 'CEL'
        model['proj']      = 'CAR'
        model['xref']      = 83.63
        model['yref']      = 22.01
        model['outcube']   = 'ctmodel_py6.fits'
        model['logfile']   = 'ctmodel_py6.log'
        model['chatter']   = 2

        # Run ctmodel tool
        model.logFileOpen()   # Make sure we get a log file
        model.execute()

        # Check result file
        self._check_result_file('ctmodel_py6.fits')

        # Setup ctmodel for stacked IRF using energy dispersion
        model = ctools.ctmodel()
        model['inobs']     = 'NONE'
        model['incube']    = 'NONE'
        model['inmodel']   = self._bkgmodel
        model['expcube']   = self._expcube
        model['psfcube']   = self._psfcube
        model['edispcube'] = self._edispcube
        model['bkgcube']   = self._bkgcube
        model['edisp']     = True
        model['rad']       = 5
        model['ra']        = 83.63
        model['dec']       = 22.01
        model['tmin']      = 0
        model['tmax']      = 1800
        model['ebinalg']   = 'LOG'
        model['emin']      = 1.0
        model['emax']      = 100.0
        model['enumbins']  = 5
        model['nxpix']     = 40
        model['nypix']     = 40
        model['binsz']     = 0.1
        model['coordsys']  = 'CEL'
        model['proj']      = 'CAR'
        model['xref']      = 83.63
        model['yref']      = 22.01
        model['outcube']   = 'ctmodel_py7.fits'
        model['logfile']   = 'ctmodel_py7.log'
        model['chatter']   = 2

        # Run ctmodel tool
        model.logFileOpen()   # Make sure we get a log file
        model.execute()

        # Check result file
        self._check_result_file('ctmodel_py7.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nx=40, ny=40, nebins=5):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Model cube file name
        nx : int, optional
            Number of pixels in X direction
        ny : int, optional
            Number of pixels in Y direction
        nebins : int, optional
            Number of energy bins
        """
        # Open result file
        cube = gammalib.GCTAEventCube(filename)

        # Check cube
        self._check_cube(cube, nx=nx, ny=ny, nebins=nebins)

        # Return
        return

    # Check model cube
    def _check_cube(self, cube, nx=40, ny=40, nebins=5):
        """
        Check model cube

        Parameters
        ----------
        cube : `~gammalib.GCTAEventCube`
            Model cube
        nx : int, optional
            Number of pixels in X direction
        ny : int, optional
            Number of pixels in Y direction
        nebins : int, optional
            Number of energy bins
        """
        # Check number of bins per axis
        self.test_value(cube.nx(), nx, 'Check number of X pixels')
        self.test_value(cube.ny(), ny, 'Check number of Y pixels')
        self.test_value(cube.ebins(), nebins, 'Check number of energy bins')

        # Return
        return
