#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comlixmap script.
#
# Copyright (C) 2021 Juergen Knoedlseder
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
# Test class for comlixmap script #
# =============================== #
class Test(test):
    """
    Test class for comlixmap script

    This test class makes unit tests for the comlixmap script by using it
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
        self._obs   = self._datadir + '/comptel/obs_binned.xml'
        self._model = self._datadir + '/comptel/models.xml'

        # Delete datastores
        os.system('rm -rf comlixmap_cmd1')
        os.system('rm -rf comlixmap_cmd2')
        os.system('rm -rf comlixmap_py1')
        os.system('rm -rf comlixmap_py2')
        os.system('rm -rf comlixmap_py3')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comlixmap')

        # Append tests
        self.append(self._test_cmd, 'Test comlixmap on command line')
        self.append(self._test_python, 'Test comlixmap from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comlixmap = self._comscript('comlixmap')

        # Setup command
        cmd = comlixmap+' inobs="'+self._obs+'" inmodel="'+self._model+'" '+ \
                          'coordsys="CEL" proj="TAN" xref="83.6331" yref="22.0145" '+ \
                          'nxpix="2" nypix="1" binsz="1.0" '+ \
                          'bkgmethod="PHINOR" srcname="Crab" '+ \
                          'outmap="comlixmap_cmd1.fits" '+ \
                          'logfile="comlixmap_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comlixmap_cmd1.fits', 2, 1)

        # Setup command
        cmd = comlixmap+' inobs="file_that_does_not_exist.xml" inmodel="'+self._model+'" '+ \
                          'coordsys="CEL" proj="TAN" xref="83.6331" yref="22.0145" '+ \
                          'nxpix="2" nypix="1" binsz="1.0" '+ \
                          'bkgmethod="PHINOR" srcname="Crab" '+ \
                          'outmap="comlixmap_cmd2.fits" '+ \
                          'logfile="comlixmap_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comlixmap)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Test PHINOR
        lix = comscripts.comlixmap()
        lix['inobs']     = self._obs
        lix['inmodel']   = self._model
        lix['srcname']   = 'Crab'
        lix['bkgmethod'] = 'PHINOR'
        lix['coordsys']  = 'CEL'
        lix['proj']      = 'TAN'
        lix['xref']      = 83.6331
        lix['yref']      = 22.0145
        lix['nxpix']     = 1
        lix['nypix']     = 2
        lix['binsz']     = 1.0
        lix['outmap']    = 'comlixmap_py1.fits'
        lix['logfile']   = 'comlixmap_py1.log'
        lix['chatter']   = 2

        # Run script and save result
        lix.logFileOpen()   # Make sure we get a log file
        lix.run()
        lix.save()

        # Check result
        self._check_result('comlixmap_py1.fits', 1, 2)

        # Test BGDLIXA
        lix = comscripts.comlixmap()
        lix['inobs']     = self._obs
        lix['inmodel']   = self._model
        lix['srcname']   = 'Crab'
        lix['bkgmethod'] = 'BGDLIXA'
        lix['coordsys']  = 'CEL'
        lix['proj']      = 'TAN'
        lix['xref']      = 83.6331
        lix['yref']      = 22.0145
        lix['nxpix']     = 1
        lix['nypix']     = 2
        lix['binsz']     = 1.0
        lix['outmap']    = 'comlixmap_py2.fits'
        lix['logfile']   = 'comlixmap_py2.log'
        lix['chatter']   = 3

        # Run script and save result
        lix.logFileOpen()   # Make sure we get a log file
        lix.execute()

        # Check result
        self._check_result('comlixmap_py1.fits', 1, 2)

        # Test BGDLIXE
        lix = comscripts.comlixmap()
        lix['inobs']     = self._obs
        lix['inmodel']   = self._model
        lix['srcname']   = 'Crab'
        lix['bkgmethod'] = 'BGDLIXE'
        lix['coordsys']  = 'CEL'
        lix['proj']      = 'TAN'
        lix['xref']      = 83.6331
        lix['yref']      = 22.0145
        lix['nxpix']     = 1
        lix['nypix']     = 2
        lix['binsz']     = 1.0
        lix['outmap']    = 'comlixmap_py3.fits'
        lix['logfile']   = 'comlixmap_py3.log'
        lix['chatter']   = 4

        # Run script and save result
        lix.logFileOpen()   # Make sure we get a log file
        lix.execute()

        # Check result
        self._check_result('comlixmap_py3.fits', 1, 2)

        # Return
        return

    # Check result file
    def _check_result(self, filename, nx, ny):
        """
        Check result file
        """
        # Load residual map file
        map = gammalib.GSkyMap(filename)

        # Check that there is one observation
        self.test_value(map.nx(), nx, 'Check for number of X pixels in sky map')
        self.test_value(map.ny(), ny, 'Check for number of Y pixels in sky map')

        # Return
        return
