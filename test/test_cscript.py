#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cspull script.
#
# Copyright (C) 2016 Juergen Knoedlseder
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
import gammalib
import ctools
from testing import test


# ================== #
# cscript_test class #
# ================== #
class cscript_test(ctools.cscript):
    """
    Test class for the cscript base class
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'cscript_test'
        self._version = '1.2.0'

        # Initialise observation container from constructor arguments
        self._obs, argv = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return

    # Check is_valid_filename() method
    def check_is_valid_filename(self, filename):
        return self._is_valid_filename(gammalib.GFilename(filename))

    # Check setup_observations() method
    def check_setup_observations(self, filename, events=True, cube=True):
        self.pars()['inobs'].value(filename)
        obs = gammalib.GObservations()
        self._setup_observations(obs, False, events, cube)
        return

    # Check setup_models() method
    def check_setup_models(self):
        obs = gammalib.GObservations()
        self._setup_models(obs)
        return

    # Check create_ebounds() method
    def check_create_ebounds(self, filename):
        self.pars()['ebinfile'].value(filename)
        return self._create_ebounds()

    # Check require_inobs() method
    def check_require_inobs(self, filename):
        self.pars()['inobs'].value(filename)
        self._require_inobs('check_require_inobs')
        return

    # Check require_inobs_nocube() method
    def check_require_inobs_nocube(self, filename):
        self.pars()['inobs'].value(filename)
        self._require_inobs_nocube('check_require_inobs_nocube')
        return

    # Check get_roi() method
    def check_get_roi(self):
        return self._get_roi()

    # Check restore_edisp() method
    def check_restore_edisp(self):
        obs   = gammalib.GObservations()
        edisp = [True, True]
        self._restore_edisp(obs, edisp)
        return

    # Check set_obs_bounds() method
    def check_set_obs_bounds(self):
        obs    = gammalib.GObservations()
        events = gammalib.GCTAEventList()
        cta    = gammalib.GCTAObservation()
        cta.events(events)
        obs.append(cta)
        self._set_obs_bounds(obs)
        return obs.copy()

    # Check get_mean_pointing() method
    def check_get_mean_pointing(self):
        obs = gammalib.GObservations()
        return self._get_mean_pointing(obs)

    # Check set_obs_response() method
    def check_set_obs_response(self, expcube, psfcube, bkgcube, edispcube, edisp):
        obs    = gammalib.GCTAObservation()
        counts = gammalib.GCTAEventCube()
        obs.events(counts)
        self.pars()['expcube'].value(expcube)
        self.pars()['psfcube'].value(psfcube)
        self.pars()['bkgcube'].value(bkgcube)
        self.pars()['edispcube'].value(edispcube)
        self.pars()['edisp'].value(edisp)
        self._set_obs_response(obs)
        return

    # Check get_current_rss() method
    def check_get_current_rss(self):
        return self._get_current_rss()


# ======================== #
# ctobservation_test class #
# ======================== #
class ctobservation_test(ctools.csobservation):
    """
    Test class for the ctobservation base class
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'ctobservation_test'
        self._version = '1.2.0'

        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(argv)

        # Return
        return

    # Check first_unbinned_observation() method
    def check_first_unbinned_observation(self):
        return self._first_unbinned_observation()

    # Check next_unbinned_observation() method
    def check_next_unbinned_observation(self):
        return self._next_unbinned_observation()


# ======================= #
# ctlikelihood_test class #
# ======================= #
class ctlikelihood_test(ctools.cslikelihood):
    """
    Test class for the ctlikelihood base class
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'ctlikelihood_test'
        self._version = '1.2.0'

        # Initialise application by calling the appropriate class constructor
        self._init_cslikelihood(argv)

        # Return
        return

    # Check evaluate() method
    def check_evaluate(self, value):
        par = gammalib.GModelPar()
        par.value(2.0)
        par.range(1.0, 3.0)
        return self._evaluate(par, value)


# ============================ #
# Test class for base classes  #
# ============================ #
class Test(test):
    """
    Test class for base classes
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Create parameter file for cscript_test application
        pars = gammalib.GApplicationPars()
        pars.append(gammalib.GApplicationPar('inobs','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('inmodel','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('ebinalg','s','h','FILE','','',''))
        pars.append(gammalib.GApplicationPar('ebinfile','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('ra','r','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('dec','r','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('rad','r','h','1.0','','',''))
        pars.append(gammalib.GApplicationPar('emin','r','h','1.0','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','10.0','','',''))
        pars.append(gammalib.GApplicationPar('expcube','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('psfcube','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('bkgcube','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('edispcube','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('edisp','b','h','yes','','',''))
        pars.append(gammalib.GApplicationPar('caldb','s','h','prod2','','',''))
        pars.append(gammalib.GApplicationPar('irf','s','h','South_0.5h','','',''))
        pars.append_standard()
        pars.save('cscript_test.par')

        # Create parameter file for ctobservation_test.par application
        pars = gammalib.GApplicationPars()
        pars.append_standard()
        pars.save('ctobservation_test.par')

        # Create parameter file for ctlikelihood_test application
        pars = gammalib.GApplicationPars()
        pars.append_standard()
        pars.save('ctlikelihood_test.par')

        # Create energy boundaries file
        ebounds = gammalib.GEbounds(1, gammalib.GEnergy(1.0,'TeV'),
                                       gammalib.GEnergy(10.0,'TeV'))
        ebounds.save('test_ebinfile_ebounds.fits[EBOUNDS]', True)
        ebounds.save('test_ebinfile_energybins.fits[ENERGYBINS]', True)
        ebounds.save('test_ebinfile_bins.fits[BINS]', True)

        # Create energy bin file
        energies = gammalib.GEnergies(2, gammalib.GEnergy(1.0,'TeV'),
                                         gammalib.GEnergy(10.0,'TeV'))
        energies.save('test_ebinfile_energies.fits[ENERGIES]', True)
        energies.save('test_ebinfile_engs.fits[ENGS]', True)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('cscript')

        # Append tests
        self.append(self._test_cscript, 'Test cscript base class')
        self.append(self._test_ctobservation, 'Test ctobservation base class')
        self.append(self._test_ctlikelihood, 'Test ctlikelihood base class')

        # Return
        return

    # Test cscript base class
    def _test_cscript(self):
        """
        Test cscript base class
        """
        # Allocate empty
        empty = cscript_test()

        # Test is_valid_filename() method
        self.test_assert(empty.check_is_valid_filename(self._model),
                         'Check model definiton XML file')
        self.test_assert(not empty.check_is_valid_filename(''),
                         'Check empty filename')
        self.test_assert(not empty.check_is_valid_filename('NONE'),
                         'Check "NONE" filename')

        # Test setup_observations() for invalid file name
        self.test_try('Check setup_observations() for "NONE"')
        try:
            empty.check_setup_observations('NONE')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test setup_observations() for acceptable event list
        self.test_try('Check setup_observations() for acceptable event list')
        try:
            empty.check_setup_observations(self._events, events=True)
            self.test_try_success()
        except ValueError:
            self.test_try_failure('Exception thrown')

        # Test setup_observations() for unacceptable event list
        self.test_try('Check setup_observations() for unacceptable event list')
        try:
            empty.check_setup_observations(self._events, events=False)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test setup_observations() for acceptable counts cube
        self.test_try('Check setup_observations() for acceptable counts cube')
        try:
            empty.check_setup_observations(self._cntcube, cube=True)
            self.test_try_success()
        except ValueError:
            self.test_try_failure('Exception thrown')

        # Test setup_observations() for unacceptable counts cube
        self.test_try('Check setup_observations() for unacceptable counts cube')
        try:
            empty.check_setup_observations(self._cntcube, cube=False)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test setup_models() for invalid model definition file name
        self.test_try('Check setup_models() for invalid model definition file name')
        try:
            empty.check_setup_models()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test create_ebounds() method for "EBOUNDS" extension
        ebounds = empty.check_create_ebounds('test_ebinfile_ebounds.fits')
        self.test_value(ebounds.size(), 1, 'Check number of energy boundaries')
        self.test_value(ebounds.emin().TeV(), 1.0, 'Check minimum energy')
        self.test_value(ebounds.emax().TeV(), 10.0, 'Check maximum energy')

        # Test create_ebounds() method for "ENERGYBINS" extension
        ebounds = empty.check_create_ebounds('test_ebinfile_energybins.fits')
        self.test_value(ebounds.size(), 1, 'Check number of energy boundaries')
        self.test_value(ebounds.emin().TeV(), 1.0, 'Check minimum energy')
        self.test_value(ebounds.emax().TeV(), 10.0, 'Check maximum energy')

        # Test create_ebounds() method for "ENERGIES" extension
        ebounds = empty.check_create_ebounds('test_ebinfile_energies.fits')
        self.test_value(ebounds.size(), 1, 'Check number of energy boundaries')
        self.test_value(ebounds.emin().TeV(), 1.0, 'Check minimum energy')
        self.test_value(ebounds.emax().TeV(), 10.0, 'Check maximum energy')

        # Test create_ebounds() method for invalid extension
        self.test_try('Check create_ebounds() method for invalid extension')
        try:
            empty.check_create_ebounds('test_ebinfile_energies.fits[INVALID]')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test create_ebounds() method for "BINS" extension
        ebounds = empty.check_create_ebounds('test_ebinfile_bins.fits[BINS]')
        self.test_value(ebounds.size(), 1, 'Check number of energy boundaries')
        self.test_value(ebounds.emin().TeV(), 1.0, 'Check minimum energy')
        self.test_value(ebounds.emax().TeV(), 10.0, 'Check maximum energy')

        # Test create_ebounds() method for "ENGS" extension
        ebounds = empty.check_create_ebounds('test_ebinfile_engs.fits[ENGS]')
        self.test_value(ebounds.size(), 1, 'Check number of energy boundaries')
        self.test_value(ebounds.emin().TeV(), 1.0, 'Check minimum energy')
        self.test_value(ebounds.emax().TeV(), 10.0, 'Check maximum energy')

        # Test require_inobs() method for invalid file name
        self.test_try('Check require_inobs() method for invalid file name')
        try:
            empty.check_require_inobs('NONE')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test require_inobs_nocube() method for valid file name
        self.test_try('Check require_inobs_nocube() method for event list')
        try:
            empty.check_require_inobs_nocube(self._events)
            self.test_try_success()
        except ValueError:
            self.test_try_failure('Exception thrown')

        # Test require_inobs_nocube() method for invalid file name
        self.test_try('Check require_inobs_nocube() method for counts cube')
        try:
            empty.check_require_inobs_nocube(self._cntcube)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test get_roi() method
        roi = empty.check_get_roi()
        self.test_assert(not roi.is_valid(), 'Check get_roi() method')

        # Test restore_edisp() method for invalid vector size
        self.test_try('Check restore_edisp() method for invalid vector size')
        try:
            empty.check_restore_edisp()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test set_obs_bounds() method
        obs = empty.check_set_obs_bounds()
        self.test_value(obs.size(), 1, 'Check number of observations')
        self.test_value(obs[0].events().ebounds().size(), 1,
                        'Check number of energy boundaries')
        self.test_value(obs[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(obs[0].events().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy')
        self.test_value(obs[0].events().roi().radius(), 1.0,
                        'Check RoI radius')

        # Test get_mean_pointing() method for empty observation container
        self.test_try('Check get_mean_pointing() method for empty observation '
                      'container')
        try:
            empty.check_get_mean_pointing()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test set_obs_response() method with 'NONE' files. If energy dispersion
        # is requested then an exception will be thrown.
        empty.check_set_obs_response('NONE', 'NONE', 'NONE', 'NONE', 'no')
        empty.check_set_obs_response('NONE', 'NONE', 'NONE', 'NONE', 'yes')
        self.test_try('Check set_obs_response() method for file names')
        try:
            empty.check_set_obs_response(self._datadir+'/crab_expcube.fits',
                                         self._datadir+'/crab_psfcube.fits',
                                         self._datadir+'/crab_bkgcube.fits',
                                         'NONE', 'yes')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()
        empty.check_set_obs_response(self._datadir+'/crab_expcube.fits',
                                     self._datadir+'/crab_psfcube.fits',
                                     self._datadir+'/crab_bkgcube.fits',
                                     self._datadir+'/crab_edispcube.fits',
                                     'yes')

        # Test get_current_rss() method
        empty.check_get_current_rss()

        # Return
        return

    # Test ctobservation base class
    def _test_ctobservation(self):
        """
        Test ctobservation base class
        """
        # Allocate empty
        empty = ctobservation_test()

        # Test first_unbinned_observation() method
        self.test_assert(empty.check_first_unbinned_observation() is None,
                         'Check first_unbinned_observation() method')

        # Test next_unbinned_observation() method
        self.test_assert(empty.check_next_unbinned_observation() is None,
                         'Check next_unbinned_observation() method')

        # Return
        return

    # Test ctlikelihood base class
    def _test_ctlikelihood(self):
        """
        Test ctlikelihood base class
        """
        # Allocate empty
        empty = ctlikelihood_test()

        # Test evaluate() method for minimum value violation
        self.test_try('Check evaluate() for minimum value violation')
        try:
            empty.check_evaluate(0.0)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test evaluate() method for maximum value violation
        self.test_try('Check evaluate() for maximum value violation')
        try:
            empty.check_evaluate(4.0)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Return
        return
