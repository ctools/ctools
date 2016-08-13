#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctpsfcube tool.
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
# Test class for ctpsfcube tool #
# ============================= #
class Test(test):
    """
    Test class for ctpsfcube tool
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
        self.name('ctpsfcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctpsfcube on command line')
        self.append(self._test_python, 'Test ctpsfcube from Python')

        # Return
        return

    # Test ctpsfcube on command line
    def _test_cmd(self):
        """
        Test ctpsfcube on the command line.
        """
        # Set tool name
        ctpsfcube = self._tool('ctpsfcube')

        # Setup ctpsfcube command
        cmd = ctpsfcube+' inobs="'+self._events+'"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctpsfcube_cmd1.fits"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=200 nypix=200 binsz=0.02'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' amax=0.3 anumbins=10'+ \
                        ' logfile="ctpsfcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctpsfcube_cmd1.fits')

        # Setup ctpsfcube command
        cmd = ctpsfcube+' inobs="event_file_that_does_not_exist.fits"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctpsfcube_cmd2.fits"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=200 nypix=200 binsz=0.02'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' amax=0.3 anumbins=10'+ \
                        ' logfile="ctpsfcube_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctpsfcube --help option
        cmd = ctpsfcube+' --help'

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

    # Test ctpsfcube from Python
    def _test_python(self):
        """
        Test ctpsfcube from Python
        """
        # Allocate ctpsfcube
        psfcube = ctools.ctpsfcube()

        # Check that empty ctpsfcube tool holds a PSF cube that has no energy
        # bins
        self._check_cube(psfcube.psfcube(), nenergies=0)

        # Check that saving does not nothing
        psfcube['outcube'] = 'ctpsfcube_py0.fits'
        psfcube['logfile'] = 'ctpsfcube_py0.log'
        psfcube.logFileOpen()
        psfcube.save()
        self.test_assert(not os.path.isfile('ctpsfcube_py0.fits'),
             'Check that no PSF cube has been created')

        # Check that clearing does not lead to an exception or segfault
        psfcube.clear()

        # Now set ctpsfcube parameters
        psfcube['inobs']    = self._events
        psfcube['incube']   = 'NONE'
        psfcube['caldb']    = self._caldb
        psfcube['irf']      = self._irf
        psfcube['ebinalg']  = 'LOG'
        psfcube['emin']     = 0.1
        psfcube['emax']     = 100
        psfcube['enumbins'] = 20
        psfcube['nxpix']    = 10
        psfcube['nypix']    = 10
        psfcube['binsz']    = 0.4
        psfcube['coordsys'] = 'CEL'
        psfcube['proj']     = 'CAR'
        psfcube['xref']     = 83.63
        psfcube['yref']     = 22.01
        psfcube['amax']     = 0.3
        psfcube['anumbins'] = 10
        psfcube['outcube']  = 'ctpsfcube_py1.fits'
        psfcube['logfile']  = 'ctpsfcube_py1.log'
        psfcube['chatter']  = 2

        # Run ctpsfcube tool
        psfcube.logFileOpen()   # Make sure we get a log file
        psfcube.run()
        psfcube.save()

        # Check result file
        self._check_result_file('ctpsfcube_py1.fits')

        # Copy ctpsfcube tool
        cpy_psfcube = psfcube.copy()

        # Check PSF cube of ctpsfcube copy
        self._check_cube(cpy_psfcube.psfcube())

        # Execute copy of ctpsfcube tool again, now with a higher chatter
        # level than before
        cpy_psfcube['emin']      = 0.2
        cpy_psfcube['emax']      = 150.0
        cpy_psfcube['outcube']   = 'ctpsfcube_py2.fits'
        cpy_psfcube['logfile']   = 'ctpsfcube_py2.log'
        cpy_psfcube['addbounds'] = True
        cpy_psfcube['chatter']   = 3
        cpy_psfcube.logFileOpen()  # Needed to get a new log file
        cpy_psfcube.execute()

        # Check result file
        self._check_result_file('ctpsfcube_py2.fits', nenergies=23)

        # Now clear copy of ctpsfcube tool
        cpy_psfcube.clear()

        # Check that the cleared copy has also cleared the PSF cube
        self._check_cube(cpy_psfcube.psfcube(), nenergies=0)

        # Get mixed observation container
        obs = self._obs_mixed()

        # Set-up ctpsfcube from observation container and input counts
        # cube
        psfcube = ctools.ctpsfcube(obs)
        psfcube['incube']    = self._cntcube
        psfcube['caldb']     = self._caldb
        psfcube['irf']       = self._irf
        #psfcube['ebinalg']   = 'LOG'
        #psfcube['emin']      = 0.1
        #psfcube['emax']      = 100
        #psfcube['enumbins']  = 20
        #psfcube['nxpix']     = 10
        #psfcube['nypix']     = 10
        #psfcube['binsz']     = 0.4
        #psfcube['coordsys']  = 'CEL'
        #psfcube['proj']      = 'CAR'
        #psfcube['xref']      = 83.63
        #psfcube['yref']      = 22.01
        psfcube['amax']      = 0.3
        psfcube['anumbins']  = 10
        psfcube['addbounds'] = True
        psfcube['outcube']   = 'ctpsfcube_py3.fits'
        psfcube['logfile']   = 'ctpsfcube_py3.log'
        psfcube['chatter']   = 4

        # Execute ctpsfcube tool
        psfcube.logFileOpen()   # Make sure we get a log file
        psfcube.execute()

        # Check result file
        self._check_result_file('ctpsfcube_py3.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nenergies=21):
        """
        Check result file

        Parameters
        ----------
        filename : str
            PSF cube file name
        nenergies : int, optional
            Number of energy bins
        """
        # Open PSF cube
        cube = gammalib.GCTACubePsf(filename)

        # Check cube
        self._check_cube(cube, nenergies=nenergies)

        # Return
        return

    # Check PSF cube
    def _check_cube(self, cube, nenergies=21):
        """
        Check PSF cube

        Parameters
        ----------
        cube : `~gammalib.GCTACubePsf`
            PSF cube
        nenergies : int, optional
            Number of energy bins
        """
        # Check dimensions
        self.test_value(len(cube.energies()), nenergies,
             'Check number of energy maps')

        # Return
        return
