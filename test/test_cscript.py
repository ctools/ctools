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
    def _check_is_valid_filename(self, filename):
        return self._is_valid_filename(gammalib.GFilename(filename))

    # Check setup_observations() method
    def _check_setup_observations(self, filename, list=True, cube=True):
        self.pars()['inobs'].value(filename)
        obs = gammalib.GObservations()
        self._setup_observations(obs, True, list, cube)
        return

    # Check setup_models() method
    def _check_setup_models(self):
        obs = gammalib.GObservations()
        self._setup_models(obs)
        return

    # Check create_ebounds() method
    def _check_create_ebounds(self, filename):
        self.pars()['ebinfile'].value(filename)
        return self._create_ebounds()

    # Check require_inobs() method
    def _check_require_inobs(self, filename):
        self.pars()['inobs'].value(filename)
        self._require_inobs('_check_require_inobs')
        return

    # Check get_roi() method
    def _check_get_roi(self):
        return self._get_roi()

    # Check restore_edisp() method
    def _check_restore_edisp(self):
        obs   = gammalib.GObservations()
        edisp = [True, True]
        self._restore_edisp(obs, edisp)
        return

    # Check set_obs_bounds() method
    def _check_set_obs_bounds(self):
        obs  = gammalib.GObservations()
        list = gammalib.GCTAEventList()
        cta  = gammalib.GCTAObservation()
        cta.events(list)
        obs.append(cta)
        self._set_obs_bounds(obs)
        return obs.copy()

    # Check get_mean_pointing() method
    def _check_get_mean_pointing(self):
        obs = gammalib.GObservations()
        return self._get_mean_pointing(obs)


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
    def _check_evaluate(self, value):
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
        ebounds.save('test_ebinfile_energies.fits[ENERGIES]', True)

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
        self.test_assert(empty._check_is_valid_filename(self._model),
                         'Check model definiton XML file')
        self.test_assert(not empty._check_is_valid_filename(''),
                         'Check empty filename')
        self.test_assert(not empty._check_is_valid_filename('NONE'),
                         'Check "NONE" filename')

        # Test setup_observations() for invalid file name
        self.test_try('Check setup_observations() for "NONE"')
        try:
            empty._check_setup_observations('NONE')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test setup_observations() for unacceptable event list
        self.test_try('Check setup_observations() for unacceptable event list')
        try:
            empty._check_setup_observations(self._events, list=False)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test setup_observations() for unacceptable counts cube
        self.test_try('Check setup_observations() for unacceptable counts cube')
        try:
            empty._check_setup_observations(self._cntcube, cube=False)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test setup_models() for invalid model definition file name
        self.test_try('Check setup_models() for invalid model definition file name')
        try:
            empty._check_setup_models()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test create_ebounds() method for "EBOUNDS" extension
        ebounds = empty._check_create_ebounds('test_ebinfile_ebounds.fits')
        self.test_value(ebounds.size(), 1, 'Check number of energy boundaries')
        self.test_value(ebounds.emin().TeV(), 1.0, 'Check minimum energy')
        self.test_value(ebounds.emax().TeV(), 10.0, 'Check maximum energy')

        # Test create_ebounds() method for "ENERGYBINS" extension
        ebounds = empty._check_create_ebounds('test_ebinfile_energybins.fits')
        self.test_value(ebounds.size(), 1, 'Check number of energy boundaries')
        self.test_value(ebounds.emin().TeV(), 1.0, 'Check minimum energy')
        self.test_value(ebounds.emax().TeV(), 10.0, 'Check maximum energy')

        # Test create_ebounds() method for invalid extension
        self.test_try('Check create_ebounds() method for invalid extension')
        try:
            empty._check_create_ebounds('test_ebinfile_energies.fits')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test require_inobs() method for invalid file name
        self.test_try('Check require_inobs() method for invalid file name')
        try:
            empty._check_require_inobs('NONE')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test get_roi() method
        roi = empty._check_get_roi()
        self.test_assert(not roi.is_valid(), 'Check get_roi() method')

        # Test restore_edisp() method for invalid vector size
        #self.test_try('Check restore_edisp() method for invalid vector size')
        #try:
        #    empty._check_restore_edisp()
        #    self.test_try_failure('Exception not thrown')
        #except ValueError:
        #    self.test_try_success()
        #
        # NOTE: This code does not work because of a missing list to
        #       integer vector conversion.

        # Test set_obs_bounds() method
        obs = empty._check_set_obs_bounds()
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
            empty._check_get_mean_pointing()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

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
        self.test_assert(empty._first_unbinned_observation() is None,
                         'Check first_unbinned_observation() method')

        # Test next_unbinned_observation() method
        self.test_assert(empty._next_unbinned_observation() is None,
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
            empty._check_evaluate(0.0)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test evaluate() method for maximum value violation
        self.test_try('Check evaluate() for maximum value violation')
        try:
            empty._check_evaluate(4.0)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Return
        return
