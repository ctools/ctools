#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comobsadd script.
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
# Test class for comobsadd script #
# =============================== #
class Test(test):
    """
    Test class for comobsadd script

    This test class makes unit tests for the comobsadd script by using it
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
        os.system('rm -rf comobsadd_cmd1')
        os.system('rm -rf comobsadd_cmd2')
        os.system('rm -rf comobsadd_py1')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comobsadd')

        # Append tests
        self.append(self._test_cmd, 'Test comobsadd on command line')
        self.append(self._test_python, 'Test comobsadd from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comobsadd = self._comscript('comobsadd')

        # Setup command
        cmd = comobsadd+' inobs="'+self._obs+'" inmodel="'+self._model+'" '+ \
                        'coordsys="CEL" proj="TAN" ra="83.6331" dec="22.0145" '+ \
                        'nchi="20" npsi="20" '+ \
                        'outfolder="comobsadd_cmd1" outobs="comobsadd_cmd1.xml" '+ \
                        'logfile="comobsadd_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comobsadd_cmd1.xml', 'comobsadd_cmd1')

        # Setup command
        cmd = comobsadd+' inobs="file_that_does_not_exist.xml" inmodel="'+self._model+'" '+ \
                        'coordsys="CEL" proj="TAN" ra="83.6331" dec="22.0145" '+ \
                        'nchi="20" npsi="20" '+ \
                        'outfolder="comobsadd_cmd2" outobs="comobsadd_cmd2.xml" '+ \
                        'logfile="comobsadd_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comobsadd)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        add = comscripts.comobsadd()
        add['inobs']     = self._obs
        add['inmodel']   = self._model
        add['coordsys']  = 'CEL'
        add['proj']      = 'TAN'
        add['ra']        = 83.6331
        add['dec']       = 22.0145
        add['nchi']      = 20
        add['npsi']      = 20
        add['outfolder'] = 'comobsadd_py1'
        add['outobs']    = 'comobsadd_py1.xml'
        add['logfile']   = 'comobsadd_py1.log'
        add['chatter']   = 2

        # Run script and save result
        add.logFileOpen()   # Make sure we get a log file
        add.run()
        add.save()

        # Check result
        self._check_result('comobsadd_py1.xml', 'comobsadd_py1')

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

        # Check DRX file
        fname = '%s/com_drx.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRX file "%s" exists' % fname)
        if exists:
            drx = gammalib.GCOMDri(fname)
            self.test_value(drx.size(), 64800, 'Check DRX size')
            self.test_value(drx.nchi(), 360, 'Check DRX Chi size')
            self.test_value(drx.npsi(), 180, 'Check DRX Psi size')
            self.test_value(drx.nphibar(), 1, 'Check DRX Phibar size')

        # Check DRG file
        fname = '%s/com_drg.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRG file "%s" exists' % fname)
        if exists:
            drg = gammalib.GCOMDri(fname)
            self.test_value(drg.size(), 10000, 'Check DRG size')
            self.test_value(drg.nchi(), 20, 'Check DRG Chi size')
            self.test_value(drg.npsi(), 20, 'Check DRG Psi size')
            self.test_value(drg.nphibar(), 25, 'Check DRG Phibar size')

        # Check DRE file
        fname = '%s/com_000750-001000keV_dre.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRE file "%s" exists' % fname)
        if exists:
            dre = gammalib.GCOMDri(fname)
            self.test_value(dre.size(), 10000, 'Check DRE size')
            self.test_value(dre.nchi(), 20, 'Check DRE Chi size')
            self.test_value(dre.npsi(), 20, 'Check DRE Psi size')
            self.test_value(dre.nphibar(), 25, 'Check DRE Phibar size')

        # Check DRB file
        fname = '%s/com_000750-001000keV_drb.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRB file "%s" exists' % fname)
        if exists:
            drb = gammalib.GCOMDri(fname)
            self.test_value(drb.size(), 10000, 'Check DRB size')
            self.test_value(drb.nchi(), 20, 'Check DRB Chi size')
            self.test_value(drb.npsi(), 20, 'Check DRB Psi size')
            self.test_value(drb.nphibar(), 25, 'Check DRB Phibar size')

        # Return
        return
