#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comlixfit script.
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
# Test class for comlixfit script #
# =============================== #
class Test(test):
    """
    Test class for comlixfit script

    This test class makes unit tests for the comlixfit script by using it
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
        os.system('rm -rf comlixfit_cmd1')
        os.system('rm -rf comlixfit_cmd2')
        os.system('rm -rf comlixfit_py1')
        os.system('rm -rf comlixfit_py2')
        os.system('rm -rf comlixfit_py3')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('comlixfit')

        # Append tests
        self.append(self._test_cmd, 'Test comlixfit on command line')
        self.append(self._test_python, 'Test comlixfit from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comlixfit = self._comscript('comlixfit')

        # Setup command
        cmd = comlixfit+' inobs="'+self._obs+'" inmodel="'+self._model+'" '+ \
                          'bkgmethod="PHINOR" outfolder="comlixfit_cmd1" '+ \
                          'outobs="comlixfit_cmd1_obs.xml" '+ \
                          'outmodel="comlixfit_cmd1_model.xml" '+ \
                          'logfile="comlixfit_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comlixfit_cmd1_obs.xml', 'comlixfit_cmd1_model.xml', 'comlixfit_cmd1')

        # Setup command
        cmd = comlixfit+' inobs="file_that_does_not_exist.xml" inmodel="'+self._model+'" '+ \
                          'bkgmethod="PHINOR" outfolder="comlixfit_cmd2" '+ \
                          'outobs="comlixfit_cmd2_obs.xml" '+ \
                          'outmodel="comlixfit_cmd2_model.xml" '+ \
                          'logfile="comlixfit_cmd2.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check --help option
        self._check_help(comlixfit)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Test PHINOR
        lix = comscripts.comlixfit()
        lix['inobs']     = self._obs
        lix['inmodel']   = self._model
        lix['bkgmethod'] = 'PHINOR'
        lix['outfolder'] = 'comlixfit_py1'
        lix['outobs']    = 'comlixfit_py1_obs.xml'
        lix['outmodel']  = 'comlixfit_py1_model.xml'
        lix['logfile']   = 'comlixfit_py1.log'
        lix['chatter']   = 2

        # Run script and save result
        lix.logFileOpen()   # Make sure we get a log file
        lix.run()
        lix.save()

        # Check result
        self._check_result('comlixfit_py1_obs.xml', 'comlixfit_py1_model.xml', 'comlixfit_py1')

        # Test BGDLIXA
        lix = comscripts.comlixfit()
        lix['inobs']     = self._obs
        lix['inmodel']   = self._model
        lix['bkgmethod'] = 'BGDLIXA'
        lix['outfolder'] = 'comlixfit_py2'
        lix['outobs']    = 'comlixfit_py2_obs.xml'
        lix['outmodel']  = 'comlixfit_py2_model.xml'
        lix['logfile']   = 'comlixfit_py2.log'
        lix['chatter']   = 3

        # Run script and save result
        lix.logFileOpen()   # Make sure we get a log file
        lix.execute()

        # Check result
        self._check_result('comlixfit_py2_obs.xml', 'comlixfit_py2_model.xml', 'comlixfit_py2')

        # Test BGDLIXE
        lix = comscripts.comlixfit()
        lix['inobs']     = self._obs
        lix['inmodel']   = self._model
        lix['bkgmethod'] = 'BGDLIXE'
        lix['outfolder'] = 'comlixfit_py3'
        lix['outobs']    = 'comlixfit_py3_obs.xml'
        lix['outmodel']  = 'comlixfit_py3_model.xml'
        lix['logfile']   = 'comlixfit_py3.log'
        lix['chatter']   = 4

        # Run script and save result
        lix.logFileOpen()   # Make sure we get a log file
        lix.execute()

        # Check result
        self._check_result('comlixfit_py3_obs.xml', 'comlixfit_py3_model.xml', 'comlixfit_py3')

        # Return
        return

    # Check result file
    def _check_result(self, obsname, modelname, foldername, nobs=4, nmodels=5):
        """
        Check result file
        """
        # Load observation definition file
        obs = gammalib.GObservations(obsname)

        # Check that there is one observation
        self.test_value(obs.size(), nobs, 'Check for %d observations' % nobs)

        # Load model definition file
        models = gammalib.GModels(modelname)

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

        # Check DRB file (0.75-1 MeV)
        drbname = '%s/vp0001_0_drb-srclix_psd0-110_000750-001000keV.fits' % (foldername)
        exists  = os.path.isfile(drbname)
        self.test_assert(exists, 'Check that file "%s" exists' % drbname)
        if exists:
            drb = gammalib.GCOMDri(drbname)
            self.test_value(drb.size(), 90000, 'Check DRB size')
            self.test_value(drb.nchi(), 60, 'Check DRB Chi size')
            self.test_value(drb.npsi(), 60, 'Check DRB Psi size')
            self.test_value(drb.nphibar(), 25, 'Check DRB Phibar size')

        # Check DRB file (1-3 MeV)
        drbname = '%s/vp0001_0_drb-srclix_psd0-110_001000-003000keV.fits' % (foldername)
        exists  = os.path.isfile(drbname)
        self.test_assert(exists, 'Check that file "%s" exists' % drbname)
        if exists:
            drb = gammalib.GCOMDri(drbname)
            self.test_value(drb.size(), 90000, 'Check DRB size')
            self.test_value(drb.nchi(), 60, 'Check DRB Chi size')
            self.test_value(drb.npsi(), 60, 'Check DRB Psi size')
            self.test_value(drb.nphibar(), 25, 'Check DRB Phibar size')

        # Check DRB file (3-10 MeV)
        drbname = '%s/vp0001_0_drb-srclix_psd0-110_003000-010000keV.fits' % (foldername)
        exists  = os.path.isfile(drbname)
        self.test_assert(exists, 'Check that file "%s" exists' % drbname)
        if exists:
            drb = gammalib.GCOMDri(drbname)
            self.test_value(drb.size(), 90000, 'Check DRB size')
            self.test_value(drb.nchi(), 60, 'Check DRB Chi size')
            self.test_value(drb.npsi(), 60, 'Check DRB Psi size')
            self.test_value(drb.nphibar(), 25, 'Check DRB Phibar size')

        # Check DRB file (10-30 MeV)
        drbname = '%s/vp0001_0_drb-srclix_psd0-110_010000-030000keV.fits' % (foldername)
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
