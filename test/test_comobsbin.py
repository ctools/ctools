#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comobsbin script.
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
# Test class for comobsbin script #
# =============================== #
class Test(test):
    """
    Test class for comobsbin script

    This test class makes unit tests for the comobsbin script by using it
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
        self._obs     = self._datadir + '/comptel/obs_selected.xml'
        self._ebounds = self._datadir + '/comptel/ebounds_std1.fits'

        # Delete datastores
        os.system('rm -rf comobsbin_cmd1')
        os.system('rm -rf comobsbin_cmd2')
        os.system('rm -rf comobsbin_py1')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comobsbin')

        # Append tests
        self.append(self._test_cmd, 'Test comobsbin on command line')
        self.append(self._test_python, 'Test comobsbin from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comobsbin = self._comscript('comobsbin')

        # Setup command
        cmd = comobsbin+' inobs="'+self._obs+'" ebinalg="FILE" '+ \
                        'ebinfile="'+self._ebounds+'" phase="NONE" '+ \
                        'coordsys="GAL" proj="TAN" '+ \
                        'outfolder="comobsbin_cmd1" outobs="comobsbin_cmd1.xml" '+ \
                        'logfile="comobsbin_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comobsbin_cmd1.xml', 'comobsbin_cmd1')

        # Setup command
        cmd = comobsbin+' inobs="file_that_does_not_exist.xml" ebinalg="FILE" '+ \
                        'ebinfile="'+self._ebounds+'" phase="NONE" '+ \
                        'coordsys="GAL" proj="TAN" '+ \
                        'outfolder="comobsbin_cmd2" outobs="comobsbin_cmd2.xml" '+ \
                        'logfile="comobsbin_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comobsbin)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        bin = comscripts.comobsbin()
        bin['inobs']     = self._obs
        bin['ebinalg']   = 'FILE'
        bin['ebinfile']  = self._ebounds
        bin['phase']     = 'NONE'
        bin['coordsys']  = 'GAL'
        bin['proj']      = 'TAN'
        bin['outfolder'] = 'comobsbin_py1'
        bin['outobs']    = 'comobsbin_py1.xml'
        bin['logfile']   = 'comobsbin_py1.log'
        bin['chatter']   = 2

        # Run script and save result
        bin.logFileOpen()   # Make sure we get a log file
        bin.run()
        bin.save()

        # Check result
        self._check_result('comobsbin_py1.xml', 'comobsbin_py1')

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
        fname = '%s/vp0001_0_drx.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRX file "%s" exists' % fname)
        if exists:
            drx = gammalib.GCOMDri(fname)
            self.test_value(drx.size(), 64800, 'Check DRX size')
            self.test_value(drx.nchi(), 360, 'Check DRX Chi size')
            self.test_value(drx.npsi(), 180, 'Check DRX Psi size')
            self.test_value(drx.nphibar(), 1, 'Check DRX Phibar size')

        # Check DRG file
        fname = '%s/vp0001_0_drg.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRG file "%s" exists' % fname)
        if exists:
            drg = gammalib.GCOMDri(fname)
            self.test_value(drg.size(), 160000, 'Check DRG size')
            self.test_value(drg.nchi(), 80, 'Check DRG Chi size')
            self.test_value(drg.npsi(), 80, 'Check DRG Psi size')
            self.test_value(drg.nphibar(), 25, 'Check DRG Phibar size')

        # Check DRE file
        fname = '%s/vp0001_0_dre_000750-001000keV.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRE file "%s" exists' % fname)
        if exists:
            dre = gammalib.GCOMDri(fname)
            self.test_value(dre.size(), 160000, 'Check DRE size')
            self.test_value(dre.nchi(), 80, 'Check DRE Chi size')
            self.test_value(dre.npsi(), 80, 'Check DRE Psi size')
            self.test_value(dre.nphibar(), 25, 'Check DRE Phibar size')

        # Check DRB file
        fname = '%s/vp0001_0_drb_000750-001000keV.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that DRB file "%s" exists' % fname)
        if exists:
            drb = gammalib.GCOMDri(fname)
            self.test_value(drb.size(), 160000, 'Check DRB size')
            self.test_value(drb.nchi(), 80, 'Check DRB Chi size')
            self.test_value(drb.npsi(), 80, 'Check DRB Psi size')
            self.test_value(drb.nphibar(), 25, 'Check DRB Phibar size')

        # Check IAQ file
        fname = '%s/iaq_000750-001000keV.fits' % foldername
        exists  = os.path.isfile(fname)
        self.test_assert(exists, 'Check that IAQ file "%s" exists' % fname)
        if exists:
            fits = gammalib.GFits(fname)
            iaq  = fits[0]
            self.test_value(iaq.npix(), 1375, 'Check IAQ size')
            self.test_value(iaq.naxis(), 2, 'Check IAQ dimension')
            self.test_value(iaq.naxes(0), 55, 'Check IAQ Phigeo size')
            self.test_value(iaq.naxes(1), 25, 'Check IAQ Phibar size')

        # Return
        return
