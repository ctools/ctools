#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctedispcube tool.
#
# Copyright (C) 2016 Maria Haupt
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


# =============================== #
# Test class for ctedispcube tool #
# =============================== #
class Test(test):
    """
    Test class for ctedispcube tool
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
        self.name('ctedispcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctedispcube on command line')
        self.append(self._test_python, 'Test ctedispcube from Python')

        # Return
        return

    # Test ctedispcube on command line
    def _test_cmd(self):
        """
        Test ctedispcube on the command line
        """
        # Set tool name
        ctedispcube = self._tool('ctedispcube')

        # Setup ctedispcube command
        cmd = ctedispcube+' inobs="'+self._events+'"'+ \
                          ' incube="NONE"'+ \
                          ' outcube="ctedispcube_cmd1.fits"'+ \
                          ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                          ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                          ' nxpix=10 nypix=10 binsz=0.4'+ \
                          ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                          ' migramax=2.0 migrabins=10'+ \
                          ' logfile="ctedispcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctedispcube_cmd1.fits')

        # Setup ctedispcube command
        cmd = ctedispcube+' inobs="events_that_do_not_exist.fits"'+ \
                          ' incube="NONE"'+ \
                          ' outcube="ctedispcube_cmd2.fits"'+ \
                          ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                          ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                          ' nxpix=10 nypix=10 binsz=0.4'+ \
                          ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                          ' migramax=2.0 migrabins=10'+ \
                          ' logfile="ctedispcube_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctedispcube --help option
        cmd = ctedispcube+' --help'

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

    # Test ctedispcube from Python
    def _test_python(self):
        """
        Test ctedispcube from Python
        """
        # Allocate ctedispcube
        edispcube = ctools.ctedispcube()

        # Check that empty ctedispcube tool holds an energy dispersion
        # cube that has no energy and migration bins
        self._check_cube(edispcube.edispcube(), nenergies=0, nmigras=0)

        # Check that saving does not nothing
        edispcube['outcube'] = 'ctedispcube_py0.fits'
        edispcube['logfile'] = 'ctedispcube_py0.log'
        edispcube.logFileOpen()
        edispcube.save()
        self.test_assert(not os.path.isfile('ctedispcube_py0.fits'),
             'Check that no energy dispersion cube has been created')

        # Check that clearing does not lead to an exception or segfault
        edispcube.clear()

        # Now set ctedispcube parameters
        edispcube['inobs']     = self._events
        edispcube['incube']    = 'NONE'
        edispcube['caldb']     = self._caldb
        edispcube['irf']       = self._irf
        edispcube['ebinalg']   = 'LOG'
        edispcube['emin']      = 0.1
        edispcube['emax']      = 100.0
        edispcube['enumbins']  = 20
        edispcube['nxpix']     = 10
        edispcube['nypix']     = 10
        edispcube['binsz']     = 0.4
        edispcube['coordsys']  = 'CEL'
        edispcube['proj']      = 'CAR'
        edispcube['xref']      = 83.63
        edispcube['yref']      = 22.01
        edispcube['migramax']  = 2.0
        edispcube['migrabins'] = 10
        edispcube['outcube']   = 'ctedispcube_py1.fits'
        edispcube['logfile']   = 'ctedispcube_py1.log'
        edispcube['chatter']   = 2

        # Run ctedispcube tool
        edispcube.logFileOpen()  # Make sure we get a log file
        edispcube.run()
        edispcube.save()

        # Check result file
        self._check_result_file('ctedispcube_py1.fits')

        # Copy ctedispcube tool
        cpy_edispcube = edispcube.copy()

        # Check energy dispersion cube of ctedispcube copy
        self._check_cube(cpy_edispcube.edispcube())

        # Execute copy of ctedispcube tool again, now with a higher chatter
        # level than before
        cpy_edispcube['emin']      = 0.2
        cpy_edispcube['emax']      = 150.0
        cpy_edispcube['outcube']   = 'ctedispcube_py2.fits'
        cpy_edispcube['logfile']   = 'ctedispcube_py2.log'
        cpy_edispcube['addbounds'] = True
        cpy_edispcube['chatter']   = 3
        cpy_edispcube.logFileOpen()  # Needed to get a new log file
        cpy_edispcube.execute()

        # Check result file
        self._check_result_file('ctedispcube_py2.fits', nenergies=23)

        # Now clear copy of ctedispcube tool
        cpy_edispcube.clear()

        # Check that the cleared copy has also cleared the energy dispersion
        # cube
        self._check_cube(cpy_edispcube.edispcube(), nenergies=0, nmigras=0)

        # Get mixed observation container
        obs = self._obs_mixed()

        # Set-up ctedispcube from observation container
        edispcube = ctools.ctedispcube(obs)
        edispcube['incube']    = 'NONE'
        edispcube['caldb']     = self._caldb
        edispcube['irf']       = self._irf
        edispcube['ebinalg']   = 'LOG'
        edispcube['emin']      = 0.1
        edispcube['emax']      = 100.0
        edispcube['enumbins']  = 20
        edispcube['nxpix']     = 10
        edispcube['nypix']     = 10
        edispcube['binsz']     = 0.4
        edispcube['coordsys']  = 'CEL'
        edispcube['proj']      = 'CAR'
        edispcube['xref']      = 83.63
        edispcube['yref']      = 22.01
        edispcube['migramax']  = 2.0
        edispcube['migrabins'] = 10
        edispcube['addbounds'] = True
        edispcube['outcube']   = 'ctedispcube_py3.fits'
        edispcube['logfile']   = 'ctedispcube_py3.log'
        edispcube['chatter']   = 4

        # Execute ctedispcube tool
        edispcube.logFileOpen()  # Make sure we get a log file
        edispcube.execute()

        # Check result file
        self._check_result_file('ctedispcube_py3.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nenergies=21):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Energy dispersion cube file name
        nenergies : int, optional
            Number of energy bins
        """
        # Open energy dispersion cube
        cube = gammalib.GCTACubeEdisp(filename)

        # Check cube
        self._check_cube(cube, nenergies=nenergies)

        # Return
        return

    # Check energy dispersion cube
    def _check_cube(self, cube, nenergies=21, nmigras=10):
        """
        Check energy dispersion cube

        Parameters
        ----------
        cube : `~gammalib.GCTACubeEdisp`
            Energy dispersion cube
        nenergies : int, optional
            Number of energy bins
        nmigras : int, optional
            Number of migration bins
        """
        # Check dimensions
        self.test_value(len(cube.energies()), nenergies,
             'Check number of energy maps')
        self.test_value(len(cube.migras()), nmigras,
             'Check number of migration bins')

        # Return
        return
