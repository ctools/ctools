#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctexpcube tool.
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
# Test class for ctexpcube tool #
# ============================= #
class Test(test):
    """
    Test class for ctexpcube tool
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
        self.name('ctexpcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctexpcube on command line')
        self.append(self._test_python, 'Test ctexpcube from Python')

        # Return
        return

    # Test ctexpcube on command line
    def _test_cmd(self):
        """
        Test ctexpcube on the command line
        """
        # Set tool name
        ctexpcube = self._tool('ctexpcube')

        # Setup ctexpcube command
        cmd = ctexpcube+' inobs="'+self._events+'"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctexpcube_cmd1.fits"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=200 nypix=200 binsz=0.02'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctexpcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctexpcube_cmd1.fits')

        # Setup ctexpcube command
        cmd = ctexpcube+' inobs="event_file_that_does_not_exist.fits"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctexpcube_cmd2.fits"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=200 nypix=200 binsz=0.02'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctexpcube_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctexpcube --help option
        cmd = ctexpcube+' --help'

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

    # Test ctexpcube from Python
    def _test_python(self):
        """
        Test ctexpcube from Python
        """
        # Allocate ctexpcube
        expcube = ctools.ctexpcube()

        # Check that empty ctexpcube tool holds an exposure cube that has
        # no energy bins
        self._check_cube(expcube.expcube(), nenergies=0)

        # Check that saving does not nothing
        expcube['outcube'] = 'ctexpcube_py0.fits'
        expcube['logfile'] = 'ctexpcube_py0.log'
        expcube.logFileOpen()
        expcube.save()
        self.test_assert(not os.path.isfile('ctexpcube_py0.fits'),
             'Check that no exposure cube has been created')

        # Check that clearing does not lead to an exception or segfault
        expcube.clear()

        # Now set ctexpcube parameters
        expcube['inobs']    = self._events
        expcube['incube']   = 'NONE'
        expcube['caldb']    = self._caldb
        expcube['irf']      = self._irf
        expcube['ebinalg']  = 'LOG'
        expcube['emin']     = 0.1
        expcube['emax']     = 100.0
        expcube['enumbins'] = 20
        expcube['nxpix']    = 200
        expcube['nypix']    = 200
        expcube['binsz']    = 0.02
        expcube['coordsys'] = 'CEL'
        expcube['proj']     = 'CAR'
        expcube['xref']     = 83.63
        expcube['yref']     = 22.01
        expcube['logfile']  = 'ctexpcube_py1.log'
        expcube['outcube']  = 'ctexpcube_py1.fits'
        expcube['chatter']  = 2

        # Run ctexpcube tool
        expcube.logFileOpen()   # Make sure we get a log file
        expcube.run()
        expcube.save()

        # Check result file
        self._check_result_file('ctexpcube_py1.fits')

        # Copy ctexpcube tool
        cpy_expcube = expcube.copy()

        # Check exposure cube of ctexpcube copy
        self._check_cube(cpy_expcube.expcube())

        # Execute copy of ctexpcube tool again, now with a higher chatter
        # level than before
        cpy_expcube['emin']      = 0.2
        cpy_expcube['emax']      = 150.0
        cpy_expcube['outcube']   = 'ctexpcube_py2.fits'
        cpy_expcube['logfile']   = 'ctexpcube_py2.log'
        cpy_expcube['publish']   = True
        cpy_expcube['addbounds'] = True
        cpy_expcube['chatter']   = 3
        cpy_expcube.logFileOpen()  # Needed to get a new log file
        cpy_expcube.execute()

        # Check result file
        self._check_result_file('ctexpcube_py2.fits', nenergies=23)

        # Now clear copy of ctexpcube tool
        cpy_expcube.clear()

        # Check that the cleared copy has also cleared the exposure cube
        self._check_cube(cpy_expcube.expcube(), nenergies=0)

        # Get mixed observation container
        obs = self._obs_mixed()

        # Set-up ctexpcube from observation container
        expcube = ctools.ctexpcube(obs)
        expcube['incube']   = 'NONE'
        expcube['caldb']    = self._caldb
        expcube['irf']      = self._irf
        expcube['ebinalg']  = 'LOG'
        expcube['emin']     = 0.1
        expcube['emax']     = 100.0
        expcube['enumbins'] = 20
        expcube['nxpix']    = 200
        expcube['nypix']    = 200
        expcube['binsz']    = 0.02
        expcube['coordsys'] = 'CEL'
        expcube['proj']     = 'CAR'
        expcube['xref']     = 83.63
        expcube['yref']     = 22.01
        expcube['logfile']  = 'ctexpcube_py3.log'
        expcube['outcube']  = 'ctexpcube_py3.fits'
        expcube['chatter']  = 4

        # Execute ctexpcube tool
        expcube.logFileOpen()  # Needed to get a new log file
        expcube.execute()

        # Check result file
        self._check_result_file('ctexpcube_py3.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nenergies=21):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Exposure cube file name
        nenergies : int, optional
            Number of energy bins
        """
        # Open exposure cube
        cube = gammalib.GCTACubeExposure(filename)

        # Check cube
        self._check_cube(cube, nenergies=nenergies)

        # Return
        return

    # Check exposure cube
    def _check_cube(self, cube, nenergies=21):
        """
        Check exposure cube

        Parameters
        ----------
        cube : `~gammalib.GCTACubeExposure`
            Exposure cube
        nenergies : int, optional
            Number of energy bins
        """
        # Check dimensions
        self.test_value(len(cube.energies()), nenergies,
             'Check number of energy maps')

        # Return
        return
