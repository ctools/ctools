#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comsrcdetect script.
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


# ================================== #
# Test class for comsrcdetect script #
# ================================== #
class Test(test):
    """
    Test class for comsrcdetect script

    This test class makes unit tests for the comsrcdetect script by using it
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
        self._map = self._datadir + '/comptel/tsmap.fits'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comsrcdetect')

        # Append tests
        self.append(self._test_cmd, 'Test comsrcdetect on command line')
        self.append(self._test_python, 'Test comsrcdetect from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comsrcdetect = self._comscript('comsrcdetect')

        # Setup command
        cmd = comsrcdetect+' inmap="'+self._map+'" threshold="20.0" '+ \
                           'outmodel="comsrcdetect_cmd1.xml" '+ \
                           'outmap="comsrcdetect_cmd1.fits" '+ \
                           'outds9file="comsrcdetect_cmd1.reg" '+ \
                           'logfile="comsrcdetect_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comsrcdetect_cmd1.xml', 'comsrcdetect_cmd1.fits', 'comsrcdetect_cmd1.reg')

        # Setup command
        cmd = comsrcdetect+' inmap="file_that_does_not_exist.xml" threshold="20.0" '+ \
                           'outmodel="comsrcdetect_cmd2.xml" '+ \
                           'outmap="comsrcdetect_cmd2.fits" '+ \
                           'outds9file="comsrcdetect_cmd2.reg" '+ \
                           'logfile="comsrcdetect_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comsrcdetect)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        cat = comscripts.comsrcdetect()
        cat['inmap']      = self._map
        cat['threshold']  = 20.0
        cat['outmodel']   = 'comsrcdetect_py1.xml'
        cat['outmap']     = 'comsrcdetect_py1.fits'
        cat['outds9file'] = 'comsrcdetect_py1.reg'
        cat['logfile']    = 'comsrcdetect_py1.log'
        cat['chatter']    = 2

        # Run script and save result
        cat.logFileOpen()   # Make sure we get a log file
        cat.run()
        cat.save()

        # Check result
        self._check_result('comsrcdetect_py1.xml', 'comsrcdetect_py1.fits', 'comsrcdetect_py1.reg')

        # Return
        return

    # Check result file
    def _check_result(self, modelname, mapname, regname, nsrc=1):
        """
        Check result file
        """
        # Load model definition file
        models = gammalib.GModels(modelname)

        # Check number of models
        self.test_value(models.size(), nsrc, 'Check for %d models' % nsrc)

        # Load output TS map file
        map = gammalib.GSkyMap(mapname)

        # Check that there is one observation
        self.test_value(map.nx(), 20, 'Check for number of X pixels in sky map')
        self.test_value(map.ny(), 20, 'Check for number of Y pixels in sky map')

        # Load regions file
        regions = gammalib.GSkyRegions(regname)

        # Check number of sky regions (point regions not defined!)
        #self.test_value(regions.size(), nsrc, 'Check for %d sky regions' % nsrc)
        self.test_value(regions.size(), 0, 'Check for %d sky regions' % 0)

        # Return
        return
