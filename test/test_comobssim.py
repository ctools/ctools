#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comobssim script.
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
# Test class for comobssim script #
# =============================== #
class Test(test):
    """
    Test class for comobssim script

    This test class makes unit tests for the comobssim script by using it
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
        os.system('rm -rf comobssim_cmd1')
        os.system('rm -rf comobssim_cmd2')
        os.system('rm -rf comobssim_py1')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comobssim')

        # Append tests
        self.append(self._test_cmd, 'Test comobssim on command line')
        self.append(self._test_python, 'Test comobssim from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comobssim = self._comscript('comobssim')

        # Setup command
        cmd = comobssim+' inobs="'+self._obs+'" inmodel="'+self._model+'" '+ \
                        'add="yes" '+ \
                        'outfolder="comobssim_cmd1" outobs="comobssim_cmd1.xml" '+ \
                        'logfile="comobssim_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comobssim_cmd1.xml', 'comobssim_cmd1')

        # Setup command
        cmd = comobssim+' inobs="file_that_does_not_exist.xml"  '+ \
                        'inmodel="'+self._model+'" add="yes" '+ \
                        'outfolder="comobssim_cmd2" outobs="comobssim_cmd2.xml" '+ \
                        'logfile="comobssim_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comobssim)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        sim = comscripts.comobssim()
        sim['inobs']     = self._obs
        sim['inmodel']   = self._model
        sim['add']       = True
        sim['outfolder'] = 'comobssim_py1'
        sim['outobs']    = 'comobssim_py1.xml'
        sim['logfile']   = 'comobssim_py1.log'
        sim['chatter']   = 2

        # Run script and save result
        sim.logFileOpen()   # Make sure we get a log file
        sim.run()
        sim.save()

        # Check result
        self._check_result('comobssim_py1.xml', 'comobssim_py1')

        # Return
        return

    # Check result file
    def _check_result(self, filename, foldername, nobs=1):
        """
        Check result file
        """
        # Load observation definition file
        obs = gammalib.GObservations(filename)

        # Check that there is one observation
        self.test_value(obs.size(), nobs, 'Check for %d observations' % nobs)

        # Check DRE file
        fname = '%s/vp0001_0_dre_sim-seed000001_psd0-110_000750-001000keV.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRE file "%s" exists' % fname)
        if exists:
            dre = gammalib.GCOMDri(fname)
            self.test_value(dre.size(), 90000, 'Check DRE size')
            self.test_value(dre.nchi(), 60, 'Check DRE Chi size')
            self.test_value(dre.npsi(), 60, 'Check DRE Psi size')
            self.test_value(dre.nphibar(), 25, 'Check DRE Phibar size')

        # Return
        return
