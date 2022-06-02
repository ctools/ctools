#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comobsmodel script.
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


# ================================= #
# Test class for comobsmodel script #
# ================================= #
class Test(test):
    """
    Test class for comobsmodel script

    This test class makes unit tests for the comobsmodel script by using it
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
        self._obs = self._datadir + '/comptel/obs_binned.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comobsmodel')

        # Append tests
        self.append(self._test_cmd, 'Test comobsmodel on command line')
        self.append(self._test_python, 'Test comobsmodel from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comobsmodel = self._comscript('comobsmodel')

        # Setup command
        cmd = comobsmodel+' inobs="'+self._obs+'" ra="83.6331" dec="22.0145" '+ \
                          'srcname="Crab" brems="NONE" ic="NONE"  iso="NONE" '+ \
                          'outmodel="comobsmodel_cmd1.xml" '+ \
                          'logfile="comobsmodel_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comobsmodel_cmd1.xml')

        # Setup command
        cmd = comobsmodel+' inobs="file_that_does_not_exist.xml" '+ \
                          'ra="83.6331" dec="22.0145" '+ \
                          'srcname="Crab" brems="NONE" ic="NONE"  iso="NONE" '+ \
                          'outmodel="comobsmodel_cmd2.xml" '+ \
                          'logfile="comobsmodel_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comobsmodel)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        model = comscripts.comobsmodel()
        model['inobs']    = self._obs
        model['ra']       = 83.6331
        model['dec']      = 22.0145
        model['srcname']  = 'Crab'
        model['brems']    = 'NONE'
        model['ic']       = 'NONE'
        model['iso']      = 'NONE'
        model['outmodel'] = 'comobsmodel_py1.xml'
        model['logfile']  = 'comobsmodel_py1.log'
        model['chatter']  = 2

        # Run script and save result
        model.logFileOpen()   # Make sure we get a log file
        model.run()
        model.save()

        # Check result
        self._check_result('comobsmodel_py1.xml')

        # Return
        return

    # Check result file
    def _check_result(self, filename, nmodels=5):
        """
        Check result file
        """
        # Load model definition file
        models = gammalib.GModels(filename)

        # Check models
        n = models.size()
        self.test_value(n, nmodels, 'Check for %d models' % nmodels)
        if n == 5:
            self.test_value(models[0].name(), 'Crab', 'Check for "Crab" in models')
            self.test_value(models[1].name(), 'Background_vp0001_0_000750-001000keV',
                'Check for "Background_vp0001_0_000750-001000keV" in models')
            self.test_value(models[2].name(), 'Background_vp0001_0_001000-003000keV',
                'Check for "Background_vp0001_0_001000-003000keV" in models')
            self.test_value(models[3].name(), 'Background_vp0001_0_003000-010000keV',
                'Check for "Background_vp0001_0_003000-010000keV" in models')
            self.test_value(models[4].name(), 'Background_vp0001_0_010000-030000keV',
                'Check for "Background_vp0001_0_010000-030000keV" in models')

        # Return
        return
