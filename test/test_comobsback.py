#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comobsback script.
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


# ================================ #
# Test class for comobsback script #
# ================================ #
class Test(test):
    """
    Test class for comobsback script

    This test class makes unit tests for the comobsback script by using it
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

        # Delete datastores
        os.system('rm -rf comobsback_cmd1')
        os.system('rm -rf comobsback_cmd2')
        os.system('rm -rf comobsback_py1')
        os.system('rm -rf comobsback_py2')
        os.system('rm -rf comobsback_py3')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comobsback')

        # Append tests
        self.append(self._test_cmd, 'Test comobsback on command line')
        self.append(self._test_python, 'Test comobsback from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comobsback = self._comscript('comobsback')

        # Setup command
        cmd = comobsback+' inobs="'+self._obs+'" inmodel="'+self._model+'" '+ \
                         'bkgmethod="PHINOR" '+ \
                         'outfolder="comobsback_cmd1" '+ \
                         'outobs="comobsback_cmd1.xml" '+ \
                         'logfile="comobsback_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comobsback_cmd1.xml', 'comobsback_cmd1', 'phinor')

        # Setup command
        cmd = comobsback+' inobs="file_that_does_not_exist.xml" '+ \
                         'inmodel="'+self._model+'" '+ \
                         'bkgmethod="PHINOR" '+ \
                         'outfolder="comobsback_cmd2" '+ \
                         'outobs="comobsback_cmd2.xml" '+ \
                         'logfile="comobsback_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comobsback)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Test PHINOR
        bgd = comscripts.comobsback()
        bgd['inobs']     = self._obs
        bgd['inmodel']   = self._model
        bgd['bkgmethod'] = 'PHINOR'
        bgd['outfolder'] = 'comobsback_py1'
        bgd['outobs']    = 'comobsback_py1.xml'
        bgd['logfile']   = 'comobsback_py1.log'
        bgd['chatter']   = 2

        # Run script and save result
        bgd.logFileOpen()   # Make sure we get a log file
        bgd.run()
        bgd.save()

        # Check result
        self._check_result('comobsback_py1.xml', 'comobsback_py1', 'phinor')

        # Test BGDLIXA
        bgd = comscripts.comobsback()
        bgd['inobs']     = self._obs
        bgd['inmodel']   = self._model
        bgd['bkgmethod'] = 'BGDLIXA'
        bgd['outfolder'] = 'comobsback_py2'
        bgd['outobs']    = 'comobsback_py2.xml'
        bgd['logfile']   = 'comobsback_py2.log'
        bgd['chatter']   = 3

        # Run script and save result
        bgd.logFileOpen()   # Make sure we get a log file
        bgd.execute()

        # Check result
        self._check_result('comobsback_py2.xml', 'comobsback_py2', 'bgdlixA-nr3-na3-ni13-ne0')

        # Test BGDLIXE
        bgd = comscripts.comobsback()
        bgd['inobs']     = self._obs
        bgd['inmodel']   = self._model
        bgd['bkgmethod'] = 'BGDLIXE'
        bgd['outfolder'] = 'comobsback_py3'
        bgd['outobs']    = 'comobsback_py3.xml'
        bgd['logfile']   = 'comobsback_py3.log'
        bgd['chatter']   = 4

        # Run script and save result
        bgd.logFileOpen()   # Make sure we get a log file
        bgd.execute()

        # Check result
        self._check_result('comobsback_py3.xml', 'comobsback_py3', 'bgdlixE-na3-ni13-ne0')

        # Return
        return

    # Check result file
    def _check_result(self, filename, foldername, bkgmethod, nobs=1):
        """
        Check result file
        """
        # Load observation definition file
        obs = gammalib.GObservations(filename)

        # Check that there is one observation
        self.test_value(obs.size(), nobs, 'Check for %d observations' % nobs)

        # Check DRB file
        drbname = '%s/vp0001_0_drb-%s_psd0-110_000750-001000keV.fits' % (foldername, bkgmethod)
        exists  = os.path.isfile(drbname)
        self.test_assert(exists, 'Check that file "%s" exists' % drbname)
        if exists:
            drb = gammalib.GCOMDri(drbname)
            self.test_value(drb.size(), 90000, 'Check DRB size')
            self.test_value(drb.nchi(), 60, 'Check DRB Chi size')
            self.test_value(drb.npsi(), 60, 'Check DRB Psi size')
            self.test_value(drb.nphibar(), 25, 'Check DRB Phibar size')

        # Return
        return
