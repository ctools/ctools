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
from testing import test


# ============================= #
# Test class for ctbkgcube tool #
# ============================= #
class Test(test):
    """
    Test class for ctbkgcube tool.
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
        # Set tool name
        ctbkgcube = self._tool('ctbkgcube')

        # Setup ctbkgcube command
        cmd = ctbkgcube+' inobs="'+self._events+'"'+ \
                        ' inmodel="'+self._model+'"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctbkgcube_cmd1.fits"'+ \
                        ' outmodel="ctbkgcube_cmd1.xml"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=10 nypix=10 binsz=0.4'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctbkgcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_files('ctbkgcube_cmd1')

        # Setup ctbkgcube command
        cmd = ctbkgcube+' inobs="event_file_that_does_not_exist.fits"'+ \
                        ' inmodel="'+self._model+'"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctbkgcube_cmd2.fits"'+ \
                        ' outmodel="ctbkgcube_cmd2.xml"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=10 nypix=10 binsz=0.4'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctbkgcube_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctbkgcube --help option
        cmd = ctbkgcube+' --help'

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

    # Test ctbkgcube from Python
    def _test_python(self):
        """
        Test ctbkgcube from Python
        """
        # Allocate ctbkgcube
        bkgcube = ctools.ctbkgcube()

        # Check that empty ctbkgcube tool holds a background cube that has
        # no energy bins
        self._check_cube(bkgcube.bkgcube(), nenergies=0)

        # Check that saving does not nothing
        bkgcube['outcube']  = 'ctbkgcube_py0.fits'
        bkgcube['outmodel'] = 'ctbkgcube_py0.xml'
        bkgcube['logfile']  = 'ctbkgcube_py0.log'
        bkgcube.logFileOpen()
        bkgcube.save()
        self.test_assert(not os.path.isfile('ctbkgcube_py0.fits'),
             'Check that no background cube has been created')

        # Check model container
        self._check_models(gammalib.GModels('ctbkgcube_py0.xml'), nmodels=0)

        # Check saving with empty file names
        bkgcube['outcube']  = ''
        bkgcube['outmodel'] = ''
        bkgcube['logfile']  = 'ctbkgcube_py1.log'
        bkgcube.logFileOpen()
        bkgcube.save()

        # Check saving with "none" model definiton name
        bkgcube['outcube']  = ''
        bkgcube['outmodel'] = 'NONE'
        bkgcube['logfile']  = 'ctbkgcube_py2.log'
        bkgcube.logFileOpen()
        bkgcube.save()

        # Check that publish method with user name does not lead to an
        # exception or segfault
        bkgcube.publish('My background cube')

        # Check that clearing does not lead to an exception or segfault
        bkgcube.clear()

        # Now set ctbkgcube parameters
        bkgcube['inobs']    = self._events
        bkgcube['inmodel']  = self._model
        bkgcube['incube']   = 'NONE'
        bkgcube['caldb']    = self._caldb
        bkgcube['irf']      = self._irf
        bkgcube['ebinalg']  = 'LOG'
        bkgcube['emin']     = 0.1
        bkgcube['emax']     = 100.0
        bkgcube['enumbins'] = 20
        bkgcube['nxpix']    = 10
        bkgcube['nypix']    = 10
        bkgcube['binsz']    = 0.4
        bkgcube['coordsys'] = 'CEL'
        bkgcube['proj']     = 'CAR'
        bkgcube['xref']     = 83.63
        bkgcube['yref']     = 22.01
        bkgcube['outcube']  = 'ctbkgcube_py3.fits'
        bkgcube['outmodel'] = 'ctbkgcube_py3.xml'
        bkgcube['logfile']  = 'ctbkgcube_py3.log'
        bkgcube['chatter']  = 2

        # Run ctbkgcube tool
        bkgcube.logFileOpen()   # Make sure we get a log file
        bkgcube.run()
        bkgcube.save()

        # Check result files
        self._check_result_files('ctbkgcube_py3')

        # Copy ctbkgcube tool
        cpy_bkgcube = bkgcube.copy()

        # Check background cube and model container of ctbkgcube copy
        self._check_cube(cpy_bkgcube.bkgcube())
        self._check_models(cpy_bkgcube.models())

        # Execute copy of ctbkgcube tool again, now with a higher chatter
        # level than before. In addition, use counts cube to define the
        # background cube
        cpy_bkgcube['incube']    = self._cntcube
        cpy_bkgcube['outcube']   = 'ctbkgcube_py4.fits'
        cpy_bkgcube['outmodel']  = 'ctbkgcube_py4.xml'
        cpy_bkgcube['logfile']   = 'ctbkgcube_py4.log'
        cpy_bkgcube['chatter']   = 3
        cpy_bkgcube['publish']   = True
        cpy_bkgcube['addbounds'] = True
        cpy_bkgcube.logFileOpen()   # Make sure we get a log file
        cpy_bkgcube.execute()

        # Check result files
        self._check_result_files('ctbkgcube_py4')

        # Now clear copy of ctbkgcube tool
        cpy_bkgcube.clear()

        # Check that the cleared copy has also cleared the background cube
        # and model container
        self._check_cube(cpy_bkgcube.bkgcube(), nenergies=0)
        self._check_models(cpy_bkgcube.models(), nmodels=0)

        # Get mixel observation container
        obs = self._obs_mixed()
        obs.models(gammalib.GModels(self._model))

        # Set-up ctbkgcube from observation container
        bkgcube = ctools.ctbkgcube(obs)
        bkgcube['incube']    = ''
        bkgcube['caldb']     = self._caldb
        bkgcube['irf']       = self._irf
        bkgcube['ebinalg']   = 'LOG'
        bkgcube['emin']      = 0.2
        bkgcube['emax']      = 150.0
        bkgcube['enumbins']  = 20
        bkgcube['nxpix']     = 10
        bkgcube['nypix']     = 10
        bkgcube['binsz']     = 0.4
        bkgcube['coordsys']  = 'CEL'
        bkgcube['proj']      = 'CAR'
        bkgcube['xref']      = 83.63
        bkgcube['yref']      = 22.01
        bkgcube['outcube']   = 'ctbkgcube_py5.fits'
        bkgcube['outmodel']  = 'ctbkgcube_py5.xml'
        bkgcube['logfile']   = 'ctbkgcube_py5.log'
        bkgcube['addbounds'] = True
        bkgcube['chatter']   = 4

        # Execute ctbkgcube tool
        bkgcube.logFileOpen()   # Make sure we get a log file
        bkgcube.execute()

        # Check result files
        self._check_result_files('ctbkgcube_py5', nenergies=23)

        # Return
        return

    # Check result files
    def _check_result_files(self, filename, nenergies=21, nmodels=2):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Test file name without extension
        nenergies : int, optional
            Number of energy bins
        """
        # Open background cube
        cube = gammalib.GCTACubeBackground(filename+'.fits')

        # Open model definition file
        models = gammalib.GModels(filename+'.xml')

        # Check cube
        self._check_cube(cube, nenergies=nenergies)

        # Check models
        self._check_models(models, nmodels=nmodels)

        # Return
        return

    # Check background cube
    def _check_cube(self, cube, nenergies=21):
        """
        Check background cube

        Parameters
        ----------
        cube : `~gammalib.GCTACubeBackground`
            Background cube
        nenergies : int, optional
            Number of energy bins
        """
        # Check dimensions
        self.test_value(len(cube.energies()), nenergies,
             'Check number of energy maps')

        # Return
        return

    # Check model container
    def _check_models(self, models, nmodels=2):
        """
        Check model container

        Parameters
        ----------
        models : `~gammalib.GModels`
            Model container
        """
        # Check number of models
        self.test_value(models.size(), nmodels, 'Check number of models')
        
        # Return
        return
