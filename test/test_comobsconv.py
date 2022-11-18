#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comobsconv script.
#
# Copyright (C) 2022 Juergen Knoedlseder
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
# Test class for comobsconv script #
# =============================== #
class Test(test):
    """
    Test class for comobsconv script

    This test class makes unit tests for the comobsconv script by using it
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
        self._model = self._datadir + '/comptel/models_fix.xml'

        # Delete datastores
        os.system('rm -rf comobsconv_cmd1')
        os.system('rm -rf comobsconv_cmd2')
        os.system('rm -rf comobsconv_py1')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comobsconv')

        # Append tests
        self.append(self._test_cmd, 'Test comobsconv on command line')
        self.append(self._test_python, 'Test comobsconv from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comobsconv = self._comscript('comobsconv')

        # Setup command
        cmd = comobsconv+' inobs="'+self._obs+'" inmodel="'+self._model+'" '+ \
                         'outfolder="comobsconv_cmd1" outobs="comobsconv_cmd1.xml" '+ \
                         'logfile="comobsconv_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comobsconv_cmd1.xml', 'comobsconv_cmd1')

        # Setup command
        cmd = comobsconv+' inobs="file_that_does_not_exist.xml" inmodel="'+self._model+'" '+ \
                         'outfolder="comobsconv_cmd2" outobs="comobsconv_cmd2.xml" '+ \
                         'logfile="comobsconv_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comobsconv)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        add = comscripts.comobsconv()
        add['inobs']     = self._obs
        add['inmodel']   = self._model
        add['outfolder'] = 'comobsconv_py1'
        add['outobs']    = 'comobsconv_py1.xml'
        add['logfile']   = 'comobsconv_py1.log'
        add['chatter']   = 2

        # Run script and save result
        add.logFileOpen()   # Make sure we get a log file
        add.run()
        add.save()

        # Check result
        self._check_result('comobsconv_py1.xml', 'comobsconv_py1')

        # Return
        return

    # Check result file
    def _check_result(self, filename, foldername, nobs=4):
        """
        Check result file
        """
        # Load observation definition file
        obs = gammalib.GObservations(filename)

        # Load models and extract Crab
        models = gammalib.GModels(self._model)
        crab   = models['Crab']

        # Check that there is one observation
        self.test_value(obs.size(), nobs, 'Check for %d observations' % nobs)

        # Check response
        self.test_value(obs[0].npred(crab), 2166.4808609, 'Check number of counts in 0.75-1 MeV bin')
        self.test_value(obs[1].npred(crab), 8775.9042976, 'Check number of counts in 1-3 MeV bin')
        self.test_value(obs[2].npred(crab), 3388.0663193, 'Check number of counts in 3-10 MeV bin')
        self.test_value(obs[3].npred(crab),  608.6723513, 'Check number of counts in 10-30 MeV bin')

        # Return
        return
