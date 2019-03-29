#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cscript base classes.
#
# Copyright (C) 2016-2018 Juergen Knoedlseder
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
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return

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

    # Check create_energies() method
    def check_create_energies(self, filename):
        self.pars()['ebinfile'].value(filename)
        return self._create_energies()

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
        pnt = gammalib.GCTAPointing()
        return self._get_roi(pnt)

    # Check get_ebounds() method
    def check_get_ebounds(self):
        return self._get_ebounds()

    # Check get_gti() method
    def check_get_gti(self):
        ref = gammalib.GTimeReference()
        return self._get_gti(ref)

    # Check get_pointing() method
    def check_get_pointing(self):
        return self._get_pointing()

    # Check get_skydir() method
    def check_get_skydir(self):
        return self._get_skydir()

    # Check set_outfile_name() method
    def check_set_outfile_name(self, filename):
        return self._set_outfile_name(filename)

    # Check is_stacked() method
    def check_is_stacked(self):
        return self._is_stacked()

    # Check is_onoff() method
    def check_is_onoff(self):
        return self._is_onoff()

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
    def check_get_mean_pointing(self, obs):
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

    # Check save_event_list() method
    def check_save_event_list(self, infile, evtname, gtiname, outfile):
        obs    = gammalib.GCTAObservation()
        events = gammalib.GCTAEventList()
        obs.events(events)
        obs.save(infile, True)
        self._save_event_list(obs, infile, evtname, gtiname, outfile)
        return

    # Check warn_too_few_energies() method
    def check_warn_too_few_energies(self, energies):
        return self._warn_too_few_energies(energies)

    # Check warn_xml_suffix() method
    def check_warn_xml_suffix(self, filename):
        return self._warn_xml_suffix(gammalib.GFilename(filename))


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
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return

    # Check first_unbinned_observation() method
    def check_first_unbinned_observation(self):
        return self._first_unbinned_observation()

    # Check next_unbinned_observation() method
    def check_next_unbinned_observation(self):
        return self._next_unbinned_observation()

    # Check set_obs_statistic() method for unbinned observation
    def check_set_obs_statistic_unbinned(self, statistic):
        obs    = gammalib.GObservations()
        events = gammalib.GCTAEventList()
        cta    = gammalib.GCTAObservation()
        cta.events(events)
        obs.append(cta)
        self.obs(obs)
        self._set_obs_statistic(statistic)
        return gammalib.toupper(self.obs()[0].statistic())

    # Check set_obs_statistic() method for binned observation
    def check_set_obs_statistic_binned(self, statistic):
        obs  = gammalib.GObservations()
        cube = gammalib.GCTAEventCube()
        cta  = gammalib.GCTAObservation()
        cta.events(cube)
        obs.append(cta)
        self.obs(obs)
        self._set_obs_statistic(statistic)
        return gammalib.toupper(self.obs()[0].statistic())

    # Check set_obs_statistic() method for On/Off observation
    def check_set_obs_statistic_onoff(self, statistic):
        obs = gammalib.GObservations()
        cta = gammalib.GCTAOnOffObservation()
        obs.append(cta)
        self.obs(obs)
        self._set_obs_statistic(statistic)
        return gammalib.toupper(self.obs()[0].statistic())

    # Check read_ogip_keywords() and write_ogip_keywords methods
    def check_ogip_keywords(self, hdu):
        self._read_ogip_keywords(hdu)
        self._write_ogip_keywords(hdu)
        self._read_ogip_keywords(hdu)
        return hdu


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
        # Initialise application by calling the appropriate class constructor
        self._init_cslikelihood(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return

    # Check evaluate() method
    def check_evaluate(self, value):
        par = gammalib.GModelPar()
        par.value(2.0)
        par.range(1.0, 3.0)
        return self._evaluate(par, value)


# =========================== #
# Test class for base classes #
# =========================== #
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
        pars.append(gammalib.GApplicationPar('tmin','t','h','2005-10-08T14:30:25','','',''))
        pars.append(gammalib.GApplicationPar('tmax','t','h','2005-10-08T14:58:26','','',''))
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

    # Setup cscript class
    def _setup_cscript(self):
        """
        Setup cscript test class
        """
        # Get test class
        cls = cscript_test()

        # Recover User parameters
        pars = cls.pars()

        # Create a copy of the User parameters
        cpy_pars = pars.copy()

        # Clear User parameters
        pars.clear()

        # Return
        return (cls, pars, cpy_pars)

    # Test cscript base class
    def _test_cscript(self):
        """
        Test cscript base class
        """
        # Allocate empty
        empty = cscript_test()

        # Test read_ahead() method
        # TODO: implement

        # Test time_reference() method
        # TODO: implement

        # Test get_observations() method
        # TODO: implement

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

        # Test create_ebounds() method for invalid extension name
        self.test_try('Check create_ebounds() method for invalid extension name')
        try:
            empty.check_create_ebounds('test_ebinfile_bins.fits')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test create_energies() method for "EBOUNDS" extension
        energies = empty.check_create_energies('test_ebinfile_ebounds.fits')
        self.test_value(energies.size(), 2, 'Check number of energies')
        self.test_value(energies[0].TeV(), 1.0, 'Check first energy')
        self.test_value(energies[1].TeV(), 10.0, 'Check second energy')

        # Test create_energies() method for "ENERGYBINS" extension
        energies = empty.check_create_energies('test_ebinfile_energybins.fits')
        self.test_value(energies.size(), 2, 'Check number of energies')
        self.test_value(energies[0].TeV(), 1.0, 'Check first energy')
        self.test_value(energies[1].TeV(), 10.0, 'Check second energy')

        # Test create_energies() method for "ENERGIES" extension
        energies = empty.check_create_energies('test_ebinfile_energies.fits')
        self.test_value(energies.size(), 2, 'Check number of energies')
        self.test_value(energies[0].TeV(), 1.0, 'Check first energy')
        self.test_value(energies[1].TeV(), 10.0, 'Check second energy')

        # Test create_energies() method for invalid extension
        self.test_try('Check create_energies() method for invalid extension')
        try:
            empty.check_create_energies('test_ebinfile_energies.fits[INVALID]')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test create_energies() method for "BINS" extension
        energies = empty.check_create_energies('test_ebinfile_bins.fits[BINS]')
        self.test_value(energies.size(), 2, 'Check number of energies')
        self.test_value(energies[0].TeV(), 1.0, 'Check first energy')
        self.test_value(energies[1].TeV(), 10.0, 'Check second energy')

        # Test create_energies() method for "ENGS" extension
        energies = empty.check_create_energies('test_ebinfile_engs.fits[ENGS]')
        self.test_value(energies.size(), 2, 'Check number of energies')
        self.test_value(energies[0].TeV(), 1.0, 'Check first energy')
        self.test_value(energies[1].TeV(), 10.0, 'Check second energy')

        # Test create_energies() method for invalid extension name
        self.test_try('Check create_energies() method for invalid extension name')
        try:
            empty.check_create_energies('test_ebinfile_bins.fits')
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test create_map() method
        # TODO: implement

        # Test create_cube() method
        # TODO: implement

        # Test create_cta_obs() method
        # TODO: implement

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

        # Test get_roi() method for the case that there is no "usepnt"
        # parameter and that there are no valid "ra" and "dec" parameters.
        # In that case an empty ROI should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('ra','r','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('dec','r','h','NONE','','',''))
        roi = cls.check_get_roi()
        self.test_assert(not roi.is_valid(), 'Check that get_roi() for absent '
                        '"usepnt" and invalid "ra" and "dec" parameters returns '
                        'an invalid ROI')
        cls.pars(cpy_pars)

        # Test get_roi() method for the case that there is no "usepnt"
        # parameter and that there is no valid "ra" parameter.
        # In that case an empty ROI should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('ra','r','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('dec','r','h','10.0','','',''))
        roi = cls.check_get_roi()
        self.test_assert(not roi.is_valid(), 'Check that get_roi() for absent '
                        '"usepnt" and an invalid "ra" parameter returns '
                        'an invalid ROI')
        cls.pars(cpy_pars)

        # Test get_roi() method for the case that there is no "usepnt"
        # parameter and that there is no valid "dec" parameter.
        # In that case an empty ROI should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('ra','r','h','10.0','','',''))
        pars.append(gammalib.GApplicationPar('dec','r','h','NONE','','',''))
        roi = cls.check_get_roi()
        self.test_assert(not roi.is_valid(), 'Check that get_roi() for absent '
                        '"usepnt" and an invalid "dec" parameter returns '
                        'an invalid ROI')
        cls.pars(cpy_pars)

        # Test get_roi() method for the case that there is no "usepnt"
        # parameter and that there are valid "ra", "dec" and "rad" parameters.
        # In that case a valid ROI should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('ra','r','h','10.0','','',''))
        pars.append(gammalib.GApplicationPar('dec','r','h','20.0','','',''))
        pars.append(gammalib.GApplicationPar('rad','r','h','1.0','','',''))
        roi = cls.check_get_roi()
        self.test_assert(roi.is_valid(), 'Check that get_roi() for absent '
                        '"usepnt" and valid "ra", "dec" and "rad" parameters '
                        'returns a valid ROI')
        self.test_value(roi.centre().dir().ra_deg(), 10.0, 'Check Right Ascension')
        self.test_value(roi.centre().dir().dec_deg(), 20.0, 'Check Declination')
        self.test_value(roi.radius(), 1.0, 'Check radius')
        cls.pars(cpy_pars)

        # Test get_roi() method for the case that there is a "usepnt"
        # parameter and that there is a valid "rad" parameter.
        # In that case a valid ROI should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('usepnt','b','h','yes','','',''))
        pars.append(gammalib.GApplicationPar('rad','r','h','1.0','','',''))
        roi = cls.check_get_roi()
        self.test_assert(roi.is_valid(), 'Check that get_roi() for present '
                        '"usepnt" and valid "rad" parameter returns a valid ROI')
        self.test_value(roi.centre().dir().ra_deg(), 0.0, 'Check Right Ascension')
        self.test_value(roi.centre().dir().dec_deg(), 0.0, 'Check Declination')
        self.test_value(roi.radius(), 1.0, 'Check radius')
        cls.pars(cpy_pars)

        # Test get_ebounds() method for the case that both "emin" and "emax"
        # parameters are invalid. In that case an empty energy boundary should
        # be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('emin','r','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','NONE','','',''))
        ebounds = cls.check_get_ebounds()
        self.test_value(ebounds.size(), 0, 'Check that get_ebounds() for '
                        'invalid "emin" and "emax" parameters returns empty '
                        'energy boundaries')
        cls.pars(cpy_pars)

        # Test get_ebounds() method for the case that the "emin" parameter is
        # invalid. In that case an empty energy boundary should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('emin','r','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','100.0','','',''))
        ebounds = cls.check_get_ebounds()
        self.test_value(ebounds.size(), 0, 'Check that get_ebounds() for '
                        'invalid "emin" parameter returns empty energy '
                        'boundaries')
        cls.pars(cpy_pars)

        # Test get_ebounds() method for the case that the "emax" parameter is
        # invalid. In that case an empty energy boundary should be returned.
        test, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('emin','r','h','1.0','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','NONE','','',''))
        ebounds = test.check_get_ebounds()
        self.test_value(ebounds.size(), 0, 'Check that get_ebounds() for '
                        'invalid "emax" parameter returns empty energy '
                        'boundaries')
        test.pars(cpy_pars)

        # Test get_ebounds() method for the case of valid "emin" and "emax"
        # parameters. In that case a single energy boundary should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('emin','r','h','1.0','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','100.0','','',''))
        ebounds = cls.check_get_ebounds()
        self.test_value(ebounds.size(), 1, 'Check that get_ebounds() for '
                        'valid "emin" and "emax" parameters returns empty one '
                        'energy boundary')
        self.test_value(ebounds.emin().TeV(), 1.0, 'Check minimum energy')
        self.test_value(ebounds.emax().TeV(), 100.0, 'Check maximum energy')
        cls.pars(cpy_pars)

        # Test get_gti() method for the case that both "tmin" and "tmax"
        # parameters are invalid. In that case an empty Good Time Interval
        # should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('tmin','t','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('tmax','t','h','NONE','','',''))
        gti = cls.check_get_gti()
        self.test_value(gti.size(), 0, 'Check that get_gti() for invalid "tmin" '
                        'and "tmax" parameters returns empty GTI')
        cls.pars(cpy_pars)

        # Test get_gti() method for the case that the "tmin" parameter is
        # invalid. In that case an empty Good Time Interval should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('tmin','t','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('tmax','t','h','2005-10-08T14:58:26','','',''))
        gti = cls.check_get_gti()
        self.test_value(gti.size(), 0, 'Check that get_gti() for an invalid '
                        '"tmin" parameter returns empty GTI')
        cls.pars(cpy_pars)

        # Test get_gti() method for the case that the "tmax" parameter is
        # invalid. In that case an empty Good Time Interval should be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('tmin','t','h','2005-10-08T14:30:25','','',''))
        pars.append(gammalib.GApplicationPar('tmax','t','h','NONE','','',''))
        gti = cls.check_get_gti()
        self.test_value(gti.size(), 0, 'Check that get_gti() for an invalid '
                        '"tmax" parameter returns empty GTI')
        cls.pars(cpy_pars)

        # Test get_gti() method for the case that both "tmin" and "tmax"
        # parameters are valid. In that case a valid Good Time Interval should
        # be returned.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('tmin','t','h','2005-10-08T14:30:25','','',''))
        pars.append(gammalib.GApplicationPar('tmax','t','h','2005-10-08T14:58:26','','',''))
        gti = cls.check_get_gti()
        self.test_value(gti.size(), 1, 'Check that get_gti() for valid "tmin" '
                        'and "tmax" parameters returns one GTI')
        self.test_value(gti.tstart().utc(), '2005-10-08T14:30:25', 'Check "tmin"')
        self.test_value(gti.tstop().utc(), '2005-10-08T14:58:26', 'Check "tmax"')
        cls.pars(cpy_pars)

        # Test get_pointing() method
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('ra','r','h','10.0','','',''))
        pars.append(gammalib.GApplicationPar('dec','r','h','20.0','','',''))
        pnt = cls.check_get_pointing()
        self.test_value(pnt.dir().ra_deg(), 10.0, 'Check Right Ascension')
        self.test_value(pnt.dir().dec_deg(), 20.0, 'Check Declination')
        cls.pars(cpy_pars)

        # Test get_skydir() method
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('ra','r','h','10.0','','',''))
        pars.append(gammalib.GApplicationPar('dec','r','h','20.0','','',''))
        dir = cls.check_get_skydir()
        self.test_value(dir.ra_deg(), 10.0, 'Check Right Ascension')
        self.test_value(dir.dec_deg(), 20.0, 'Check Declination')
        cls.pars(cpy_pars)

        # Test set_outfile_name() method for "fix" filename, which should
        # return the prefixed basename.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('prefix','s','h','pre_','','',''))
        filename = cls.check_set_outfile_name('fix')
        self.test_value(filename, 'pre_fix', 'Check "pre_fix"')
        cls.pars(cpy_pars)

        # Test set_outfile_name() method for "/path/to/fix" filename, which
        # should return the prefixed basename.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('prefix','s','h','pre_','','',''))
        filename = cls.check_set_outfile_name('/path/to/fix')
        self.test_value(filename, 'pre_fix', 'Check "pre_fix"')
        cls.pars(cpy_pars)

        # Test set_outfile_name() method for "/path/to/fix.gz" filename, which
        # should return the prefixed basename without the ".gz"
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('prefix','s','h','pre_','','',''))
        filename = cls.check_set_outfile_name('/path/to/fix.gz')
        self.test_value(filename, 'pre_fix', 'Check "pre_fix"')
        cls.pars(cpy_pars)

        # Test set_outfile_name() method for empty filename, which should
        # return an empty filename.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('prefix','s','h','pre_','','',''))
        filename = cls.check_set_outfile_name('')
        self.test_value(filename, '', 'Check on empty string')
        cls.pars(cpy_pars)

        # Test is_stacked() method for the case that "enumbins" is zero. In
        # that case the is_stacked() method should return False.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('emin','r','h','1.0','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','100.0','','',''))
        pars.append(gammalib.GApplicationPar('enumbins','i','h','0','','',''))
        self.test_assert(not cls.check_is_stacked(), 'Check if is_stacked() '
                         'returns False if there are no energy bins')
        cls.pars(cpy_pars)

        # Test is_stacked() method for the case that "enumbins" is larger
        # than zero. In that case the is_stacked() method should return True
        # and query a bunch of parameters.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('emin','r','h','1.0','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','100.0','','',''))
        pars.append(gammalib.GApplicationPar('enumbins','i','h','10','','',''))
        pars.append(gammalib.GApplicationPar('coordsys','s','h','CEL','','',''))
        pars.append(gammalib.GApplicationPar('proj','s','h','TAN','','',''))
        pars.append(gammalib.GApplicationPar('xref','r','h','11.1','','',''))
        pars.append(gammalib.GApplicationPar('yref','r','h','12.1','','',''))
        pars.append(gammalib.GApplicationPar('nxpix','i','h','100','','',''))
        pars.append(gammalib.GApplicationPar('nypix','i','h','100','','',''))
        pars.append(gammalib.GApplicationPar('binsz','r','h','0.02','','',''))
        self.test_assert(cls.check_is_stacked(), 'Check if is_stacked() '
                         'returns True if there are energy bins')
        cls.pars(cpy_pars)

        # Test is_onoff() method for the case that method='NONE', i.e. the
        # method is not 'ONOFF'. In that case the is_onoff() method should
        # return False.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('method','s','h','NONE','','',''))
        self.test_assert(not cls.check_is_onoff(),
             'Check if is_onoff() returns False if method is not "ONOFF"')
        cls.pars(cpy_pars)

        # Test is_onoff() method for the case that method='ONOFF',
        # srcshape='CIRCLE' and bkgmethod='REFLECTED'. In that case the
        # is_onoff() method should return True and a bunch of application
        # parameters are queried.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('method','s','h','ONOFF','','',''))
        pars.append(gammalib.GApplicationPar('inexclusion','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('emin','r','h','0.1','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','100.0','','',''))
        pars.append(gammalib.GApplicationPar('enumbins','i','h','10','','',''))
        pars.append(gammalib.GApplicationPar('coordsys','s','h','CEL','','',''))
        pars.append(gammalib.GApplicationPar('xref','r','h','11.1','','',''))
        pars.append(gammalib.GApplicationPar('yref','r','h','12.1','','',''))
        pars.append(gammalib.GApplicationPar('srcshape','s','h','CIRCLE','','',''))
        pars.append(gammalib.GApplicationPar('rad','r','h','1.1','','',''))
        pars.append(gammalib.GApplicationPar('bkgmethod','s','h','REFLECTED','','',''))
        pars.append(gammalib.GApplicationPar('bkgregmin','i','h','3','','',''))
        pars.append(gammalib.GApplicationPar('maxoffset','r','h','0.1','','',''))
        pars.append(gammalib.GApplicationPar('etruemin','r','h','0.1','','',''))
        pars.append(gammalib.GApplicationPar('etruemax','r','h','100.0','','',''))
        pars.append(gammalib.GApplicationPar('etruebins','i','h','10','','',''))
        self.test_assert(cls.check_is_onoff(),
             'Check if is_onoff() returns True if method is "ONOFF"')
        cls.pars(cpy_pars)

        # Test is_onoff() method for the case that method='ONOFF',
        # srcshape='CIRCLE' and bkgmethod='REFLECTED'. In that case the
        # is_onoff() method should return True and a bunch of application
        # parameters are queried, except of rad and bkgregmin.
        cls, pars, cpy_pars = self._setup_cscript()
        pars.append(gammalib.GApplicationPar('method','s','h','ONOFF','','',''))
        pars.append(gammalib.GApplicationPar('inexclusion','f','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('emin','r','h','0.1','','',''))
        pars.append(gammalib.GApplicationPar('emax','r','h','100.0','','',''))
        pars.append(gammalib.GApplicationPar('enumbins','i','h','10','','',''))
        pars.append(gammalib.GApplicationPar('coordsys','s','h','CEL','','',''))
        pars.append(gammalib.GApplicationPar('xref','r','h','11.1','','',''))
        pars.append(gammalib.GApplicationPar('yref','r','h','12.1','','',''))
        pars.append(gammalib.GApplicationPar('srcshape','s','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('bkgmethod','s','h','NONE','','',''))
        pars.append(gammalib.GApplicationPar('maxoffset','r','h','0.1','','',''))
        pars.append(gammalib.GApplicationPar('etruemin','r','h','0.1','','',''))
        pars.append(gammalib.GApplicationPar('etruemax','r','h','100.0','','',''))
        pars.append(gammalib.GApplicationPar('etruebins','i','h','10','','',''))
        self.test_assert(cls.check_is_onoff(),
             'Check if is_onoff() returns True if method is "ONOFF"')
        cls.pars(cpy_pars)

        # Test set_response() method
        # TODO: implement

        # Test set_edisp() method
        # TODO: implement

        # Test restore_edisp() method for invalid vector size
        self.test_try('Check restore_edisp() method for invalid vector size')
        try:
            empty.check_restore_edisp()
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
            obs = gammalib.GObservations()
            empty.check_get_mean_pointing(obs)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test get_mean_pointing() method for observation container with LAT obs
        self.test_try('Check get_mean_pointing() method for empty observation '
                      'container')
        try:
            lat = gammalib.GLATObservation()
            obs = gammalib.GObservations()
            obs.append(lat)
            empty.check_get_mean_pointing(obs)
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test get_current_rss() method
        empty.check_get_current_rss()

        # Test get_obs_header() method
        # TODO: implement

        # Test insert_energy_boundaries() method
        # TODO: implement

        # Test cube_layer_usage() method
        # TODO: implement

        # Check save_event_list() method
        empty.check_save_event_list('test_event_list_input1.fits', 'EVENTS',
                                    'GTI', 'test_event_list_output1.fits')
        self._test_cscript_save_event_list('test_event_list_output1.fits',
                                           ['Primary','EVENTS','GTI'])
        empty.check_save_event_list('test_event_list_input2.fits', 'EVENTS',
                                    'GTI', 'test_event_list_output2.fits[EVT]')
        self._test_cscript_save_event_list('test_event_list_output2.fits',
                                           ['Primary','EVT','GTI'])
        empty.check_save_event_list('test_event_list_input3.fits', 'EVENTS',
                                    'GTI', 'test_event_list_output3.fits[EVT;TIME]')
        self._test_cscript_save_event_list('test_event_list_output3.fits',
                                           ['Primary','EVT','TIME'])
        empty.check_save_event_list('test_event_list_input4.fits', 'EVENTS',
                                    'GTI', 'test_event_list_output4.fits[EVT;]')
        self._test_cscript_save_event_list('test_event_list_output4.fits',
                                           ['Primary','EVT','GTI'])
        empty.check_save_event_list('test_event_list_input5.fits', 'EVENTS',
                                    'GTI', 'test_event_list_output5.fits[;TIME]')
        self._test_cscript_save_event_list('test_event_list_output5.fits',
                                           ['Primary','EVENTS','TIME'])
        empty.check_save_event_list('test_event_list_input6.fits', 'XEVENTS',
                                    'XGTI', 'test_event_list_output6.fits[EVT;TIME]')
        self._test_cscript_save_event_list('test_event_list_output6.fits',
                                           ['Primary','EVT','TIME','EVENTS','GTI'])

        # Test get_gtiname() method
        # TODO: implement

        # Test warn_too_few_energies() method
        energies_empty = gammalib.GEnergies()
        energies_2     = gammalib.GEnergies(2, gammalib.GEnergy(1.0,'TeV'),
                                               gammalib.GEnergy(10.0,'TeV'))
        energies_100   = gammalib.GEnergies(100, gammalib.GEnergy(1.0,'TeV'),
                                                 gammalib.GEnergy(10.0,'TeV'))
        self.test_value(len(empty.check_warn_too_few_energies(energies_empty)), 0,
                        'Check length of warning string for too few energies')
        self.test_value(len(empty.check_warn_too_few_energies(energies_2)), 253,
                        'Check length of warning string for too few energies')
        self.test_value(len(empty.check_warn_too_few_energies(energies_100)), 0,
                        'Check length of warning string for too few energies')

        # Test warn_xml_suffix() method
        self.test_value(len(empty.check_warn_xml_suffix('test.xml')), 0,
                        'Check xml suffix warning for "test.xml" file')
        self.test_value(len(empty.check_warn_xml_suffix('test.fits')), 235,
                        'Check xml suffix warning for "test.fits" file')

        # Return
        return

    # Test check_save_event_list() method result
    def _test_cscript_save_event_list(self, filename, hdunames):
        """
        Test check_save_event_list() method result
        """
        # Open FITS file
        fits = gammalib.GFits(filename)

        # Check number of HDUs
        self.test_value(fits.size(), len(hdunames), 'save_event_list(): '
                        'Check number of FITS HDUs in file %s' % filename)

        # Check HDU extension names
        for i, hdu in enumerate(fits):
            self.test_value(hdu.extname(), hdunames[i], 'save_event_list(): '
                            'Check FITS HDU extension name in file %s' % filename)

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

        # Test set_obs_statistic() method for unbinned observation
        self.test_value(empty.check_set_obs_statistic_unbinned('CSTAT'),
                        'CSTAT',
                        'Check set_obs_statistic("CSTAT") for unbinned')
        self.test_value(empty.check_set_obs_statistic_unbinned('WSTAT'),
                        'CSTAT',
                        'Check set_obs_statistic("WSTAT") for unbinned')
        self.test_value(empty.check_set_obs_statistic_unbinned('CHI2'),
                        'CSTAT',
                        'Check set_obs_statistic("CHI2") for unbinned')
        self.test_value(empty.check_set_obs_statistic_unbinned('DEFAULT'),
                        'CSTAT',
                        'Check set_obs_statistic("DEFAULT") for unbinned')

        # Test set_obs_statistic() method for binned observation
        self.test_value(empty.check_set_obs_statistic_binned('CSTAT'),
                        'CSTAT',
                        'Check set_obs_statistic("CSTAT") for binned')
        self.test_value(empty.check_set_obs_statistic_binned('WSTAT'),
                        'CSTAT',
                        'Check set_obs_statistic("WSTAT") for binned')
        self.test_value(empty.check_set_obs_statistic_binned('CHI2'),
                        'CHI2',
                        'Check set_obs_statistic("CHI2") for binned')
        self.test_value(empty.check_set_obs_statistic_binned('DEFAULT'),
                        'CSTAT',
                        'Check set_obs_statistic("DEFAULT") for binned')

        # Test set_obs_statistic() method for On/Off observation
        self.test_value(empty.check_set_obs_statistic_onoff('CSTAT'),
                        'CSTAT',
                        'Check set_obs_statistic("CSTAT") for On/Off')
        self.test_value(empty.check_set_obs_statistic_onoff('WSTAT'),
                        'WSTAT',
                        'Check set_obs_statistic("WSTAT") for On/Off')
        self.test_value(empty.check_set_obs_statistic_onoff('CHI2'),
                        'CSTAT',
                        'Check set_obs_statistic("CHI2") for On/Off')
        self.test_value(empty.check_set_obs_statistic_onoff('DEFAULT'),
                        'CSTAT',
                        'Check set_obs_statistic("DEFAULT") for On/Off')

        # Test read_ogip_keywords() method for empty header
        hdu = empty.check_ogip_keywords(gammalib.GFitsBinTable())
        self.test_value(hdu.string('TELESCOP'), '', 'Check "TELESCOP" keyword')
        self.test_value(hdu.string('DATE-OBS'), '2010-01-01', 'Check "DATE-OBS" keyword')
        self.test_value(hdu.string('TIME-OBS'), '00:00:00', 'Check "TIME-OBS" keyword')
        self.test_value(hdu.string('DATE-END'), '2010-01-01', 'Check "DATE-END" keyword')
        self.test_value(hdu.string('TIME-END'), '00:00:00', 'Check "TIME-END" keyword')
        self.test_value(hdu.real('TELAPSE'), 0.0, 'Check "TELAPSE" keyword')
        self.test_value(hdu.real('ONTIME'), 0.0, 'Check "ONTIME" keyword')
        self.test_value(hdu.real('LIVETIME'), 0.0, 'Check "LIVETIME" keyword')
        self.test_value(hdu.real('EXPOSURE'), 0.0, 'Check "EXPOSURE" keyword')
        self.test_value(hdu.real('DEADC'), 1.0, 'Check "DEADC" keyword')

        # Test read_ogip_keywords() method for filled header
        table = gammalib.GFitsBinTable()
        table.card('TELESCOP', 'CTA', 'Comment')
        table.card('DATE-OBS', '2018-01-01', 'Comment')
        table.card('TIME-OBS', '01:02:03', 'Comment')
        table.card('DATE-END', '2018-02-05', 'Comment')
        table.card('TIME-END', '21:12:05', 'Comment')
        table.card('TELAPSE', 1234.0, 'Comment')
        table.card('ONTIME', 2345.0, 'Comment')
        table.card('LIVETIME', 2210.0, 'Comment')
        table.card('EXPOSURE', 3001.0, 'Comment')

        hdu = empty.check_ogip_keywords(table)
        self.test_value(hdu.string('TELESCOP'), 'CTA', 'Check "TELESCOP" keyword')
        self.test_value(hdu.string('DATE-OBS'), '2018-01-01', 'Check "DATE-OBS" keyword')
        self.test_value(hdu.string('TIME-OBS'), '01:02:03', 'Check "TIME-OBS" keyword')
        self.test_value(hdu.string('DATE-END'), '2018-02-05', 'Check "DATE-END" keyword')
        self.test_value(hdu.string('TIME-END'), '21:12:05', 'Check "TIME-END" keyword')
        self.test_value(hdu.real('TELAPSE'), 1234.0, 'Check "TELAPSE" keyword')
        self.test_value(hdu.real('ONTIME'), 2345.0, 'Check "ONTIME" keyword')
        self.test_value(hdu.real('LIVETIME'), 2210.0, 'Check "LIVETIME" keyword')
        self.test_value(hdu.real('EXPOSURE'), 3001.0, 'Check "EXPOSURE" keyword')
        self.test_value(hdu.real('DEADC'), 2210.0/2345.0, 'Check "DEADC" keyword')

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
