#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comobsres script.
#
# Copyright (C) 2021-2022 Juergen Knoedlseder
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
import comscripts
from testing import test


# =============================== #
# Test class for comobsres script #
# =============================== #
class Test(test):
    """
    Test class for comobsres script

    This test class makes unit tests for the comobsres script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test datasets and parameters
        self._obs   = self._datadir + '/comptel/obs_binned1.xml'
        self._model = self._datadir + '/comptel/models.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comobsres')

        # Append tests
        self.append(self._test_cmd, 'Test comobsres on command line')
        self.append(self._test_python, 'Test comobsres from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comobsres = self._comscript('comobsres')

        # Setup command
        cmd = comobsres+' inobs="'+self._obs+'" inmodel="'+self._model+'" '+ \
                        'algorithm="SIGNIFICANCE" margin="0.0" binsz="5.0" '+ \
                        'coordsys="GAL" proj="TAN" '+ \
                        'outfolder="comobsres_cmd1" outmap="comobsres_cmd1.fits" '+ \
                        'logfile="comobsres_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comobsres_cmd1.fits')

        # Setup command
        cmd = comobsres+' inobs="file_that_does_not_exist.xml" inmodel="'+self._model+'" '+ \
                        'algorithm="SIGNIFICANCE" margin="0.0" binsz="5.0" '+ \
                        'coordsys="GAL" proj="TAN" '+ \
                        'outfolder="comobsres_cmd2" outmap="comobsres_cmd2.fits" '+ \
                        'logfile="comobsres_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comobsres)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        add = comscripts.comobsres()
        add['inobs']     = self._obs
        add['inmodel']   = self._model
        add['algorithm'] = 'SIGNIFICANCE'
        add['coordsys']  = 'GAL'
        add['proj']      = 'TAN'
        add['margin']    = 0.0
        add['binsz']     = 5.0
        add['outfolder'] = 'comobsres_py1'
        add['outmap']    = 'comobsres_py1.fits'
        add['logfile']   = 'comobsres_py1.log'
        add['chatter']   = 2

        # Run script and save result
        add.logFileOpen()   # Make sure we get a log file
        add.run()
        add.save()

        # Check result
        self._check_result('comobsres_py1.fits')

        # Return
        return

    # Check result file
    def _check_result(self, filename):
        """
        Check result file
        """
        # Load residual map file
        map = gammalib.GSkyMap(filename)

        # Check that there is one observation
        self.test_value(map.npix(), 121, 'Check for number of pixels in sky map')

        # Return
        return
