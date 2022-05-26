#! /usr/bin/env python
# ==========================================================================
# Computes the array sensitivity using the Test Statistic for a test source
#
# Copyright (C) 2011-2022 Juergen Knoedlseder
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
import sys
import math
import gammalib
import ctools
from cscripts import obsutils
from cscripts import modutils
from cscripts import mputils


# ============ #
# cssens class #
# ============ #
class cssens(ctools.csobservation):
    """
    Computes the CTA sensitivity

    This class computes the CTA sensitivity for a number of energy bins using
    ctlike. Spectra are fitted in narrow energy bins to simulated data,
    and the flux level is determined that leads to a particular significance.

    The significance is determined using the Test Statistic, defined as twice
    the likelihood difference between fitting with and without the test source.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise class members
        self._ebounds     = gammalib.GEbounds()
        self._obs_ebounds = []
        self._fits        = None
        self._srcname     = ''
        self._ra          = None
        self._dec         = None
        self._log_clients = False
        self._models      = gammalib.GModels()
        self._seed        = 1
        self._nthreads    = 0

        # Return
        return

    # State methods por pickling
    def __getstate__(self):
        """
        Extend ctools.csobservation getstate method to include some members

        Returns
        -------
        state : dict
            Pickled instance
        """
        # Set pickled dictionary
        state = {'base'         : ctools.csobservation.__getstate__(self),
                 'ebounds'      : self._ebounds,
                 'obs_ebounds'  : self._obs_ebounds,
                 'fits'         : self._fits,
                 'srcname'      : self._srcname,
                 'ra'           : self._ra,
                 'dec'          : self._dec,
                 'log_clients'  : self._log_clients,
                 'models'       : self._models,
                 'seed'         : self._seed,
                 'nthreads'     : self._nthreads}

        # Return pickled dictionary
        return state

    def __setstate__(self, state):
        """
        Extend ctools.csobservation setstate method to include some members

        Parameters
        ----------
        state : dict
            Pickled instance
        """
        ctools.csobservation.__setstate__(self, state['base'])
        self._ebounds     = state['ebounds']
        self._obs_ebounds = state['obs_ebounds']
        self._fits        = state['fits']
        self._srcname     = state['srcname']
        self._ra          = state['ra']
        self._dec         = state['dec']
        self._log_clients = state['log_clients']
        self._models      = state['models']
        self._seed        = state['seed']
        self._nthreads    = state['nthreads']

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get user parameters from parfile
        """
        # Set observation if not done before
        if self.obs().size() == 0:
            self.obs(self._set_obs(self['emin'].real(), self['emax'].real()))

        # Set observation statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

        # Set models if we have none
        if self.obs().models().size() == 0:
            self.obs().models(self['inmodel'].filename())

        # Set source name and position
        self._set_source()

        # Read further parameters
        emin = self['emin'].real()
        emax = self['emax'].real()
        bins = self['bins'].integer()

        # Query parameters for binned if requested
        enumbins = self['enumbins'].integer()
        if not enumbins == 0:
            self['npix'].integer()
            self['binsz'].real()

        # Query input parameters
        self['sigma'].real()
        self['max_iter'].integer()
        self['type'].string()
        self['edisp'].boolean()
        self['debug'].boolean()
        self['mincounts'].integer()

        # Read seed
        self._seed = self['seed'].integer()

        # Derive some parameters
        self._ebounds = gammalib.GEbounds(bins,
                                          gammalib.GEnergy(emin, 'TeV'),
                                          gammalib.GEnergy(emax, 'TeV'))

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Set number of processes for multiprocessing
        self._nthreads = mputils.nthreads(self)

        # Return
        return

    def _median(self, array):
        """
        Compute median value for an array

        Parameters
        ----------
        array : list of floats
            Array

        Returns
        -------
        median : float
            Median value of array
        """
        # Initialise median
        median = 0.0

        # Get number of elements in array
        n = len(array)

        # Continue only if there are elements in the array
        if n > 0:

            # Sort the array
            array.sort()

            # Get median
            if n % 2 != 0:
                median = array[int(n/2)]
            else:
                median = (array[int((n-1)/2)] + array[int(n/2)]) / 2.0

        # Return median
        return median

    def _set_source(self):
        """
        Set source name and position
        """
        # Set source name
        self._srcname = self['srcname'].string()

        # Set source position. If the test source has no position then set the
        # source position to (RA,Dec)=(0,0)
        source = self.obs().models()[self._srcname]
        if source.has_par('RA') and source.has_par('DEC'):
            self._ra  = source['RA'].value()
            self._dec = source['DEC'].value()
        elif source.has_par('GLON') and source.has_par('GLAT'):
            glon      = source['GLON'].value()
            glat      = source['GLAT'].value()
            srcdir    = gammalib.GSkyDir()
            srcdir.lb_deg(glon, glat)
            self._ra  = srcdir.ra_deg()
            self._dec = srcdir.dec_deg()
        else:
            self._ra  = 0.0
            self._dec = 0.0

        # Return
        return

    def _set_obs(self, emin, emax):
        """
        Set an observation container

        Parameters
        ----------
        emin : float
            Minimum energy (TeV)
        emax : float
            Maximum energy (TeV)

        Returns
        -------
        obs : `~gammalib.GObservations`
            Observation container
        """
        # If an observation was provided on input then load it from XML file
        if self['inobs'].is_valid():
            obs = self._get_observations()

        # ... otherwise allocate a single observation using the test source
        # position as pointing direction, optionally offset by a certain
        # amount
        else:

            # Load models
            models = gammalib.GModels(self['inmodel'].filename())

            # Get test source
            source = models[self['srcname'].string()]

            # Set pointing direction to test source position. If test source
            # has no position then set the pointing to (RA,Dec)=(0,0)
            pntdir = gammalib.GSkyDir()
            if source.has_par('RA') and source.has_par('DEC'):
                pntdir.radec_deg(source['RA'].value(), source['DEC'].value())
            elif source.has_par('GLON') and source.has_par('GLAT'):
                pntdir.lb_deg(source['GLON'].value(), source['GLAT'].value())
            else:
                pntdir.radec_deg(0.0, 0.0)

            # Read other relevant user parameters
            instrument = self['instrument'].string()
            caldb      = self['caldb'].string()
            irf        = self['irf'].string()
            deadc      = self['deadc'].real()
            duration   = self['duration'].real()
            rad        = self['rad'].real()
            offset     = self['offset'].real()

            # Add offset to pointing direction
            pntdir.rotate_deg(0.0, offset)

            # Allocate observation container
            obs = gammalib.GObservations()

            # Create CTA observation
            run = obsutils.set_obs(pntdir, instrument=instrument, caldb=caldb, irf=irf,
                                   duration=duration, deadc=deadc,
                                   emin=emin, emax=emax, rad=rad)

            # Append observation to container
            obs.append(run)

        # Return observation container
        return obs

    def _set_obs_ebounds(self, emin, emax):
        """
        Set energy boundaries for all observations in container

        Parameters
        ----------
        emin : `~gammalib.GEnergy`
            Minimum energy
        emax : `~gammalib.GEnergy`
            Maximum energy
        """
        # Loop over all observations in container
        for i, obs in enumerate(self.obs()):

            # Get energy boundaries of the observation
            obs_ebounds = self._obs_ebounds[i]

            # Get minimum and maximum energy of the observation
            obs_emin = obs_ebounds.emin()
            obs_emax = obs_ebounds.emax()

            # If [emin,emax] is fully contained in the observation energy range
            # the use [emin,emax] as energy boundaries
            if obs_emin <= emin and obs_emax >= emax:
                ebounds = gammalib.GEbounds(emin, emax)

            # ... otherwise, if [emin,emax] is completely outside the
            # observation energy range then set the energy boundaries to the
            # zero-width interval [0,0]
            elif emax < obs_emin or emin > obs_emax:
                e0      = gammalib.GEnergy(0.0, 'TeV')
                ebounds = gammalib.GEbounds(e0, e0)

            # ... otherwise, if [emin,emax] overlaps partially with the
            # observation energy range then set the energy boundaries to the
            # overlapping part
            else:
                # Set overlapping energy range
                set_emin = max(emin, obs_emin)
                set_emax = min(emax, obs_emax)

                # Set energy boundaries
                ebounds = gammalib.GEbounds(set_emin, set_emax)

            # Set the energy boundaries as the boundaries of the observation
            obs.events().ebounds(ebounds)

        # Return
        return

    def _get_crab_flux(self, emin, emax):
        """
        Return Crab photon flux in a given energy interval

        Parameters
        ----------
        emin : `~gammalib.GEnergy`
            Minimum energy
        emax : `~gammalib.GEnergy`
            Maximum energy

        Returns
        -------
        flux : float
            Crab photon flux in specified energy interval (ph/cm2/s)
        """
        # Set Crab TeV spectral model based on a power law
        crab = gammalib.GModelSpectralPlaw(5.7e-16, -2.48,
                                           gammalib.GEnergy(0.3, 'TeV'))

        # Determine photon flux
        flux = crab.flux(emin, emax)

        # Return photon flux
        return flux

    def _get_sensitivity(self, emin, emax, test_model):
        """
        Determine sensitivity for given observations

        Parameters
        ----------
        emin : `~gammalib.GEnergy`
            Minimum energy for fitting and flux computation
        emax : `~gammalib.GEnergy`
            Maximum energy for fitting and flux computation
        test_model : `~gammalib.GModels`
            Test source model

        Returns
        -------
        result : dict
            Result dictionary
        """
        # Set TeV->erg conversion factor
        tev2erg = 1.6021764

        # Set parameters
        ts_thres = self['sigma'].real() * self['sigma'].real()
        max_iter = self['max_iter'].integer()
        enumbins = self['enumbins'].integer()
        if not enumbins == 0:
            npix  = self['npix'].integer()
            binsz = self['binsz'].real()
        else:
            npix  = 200
            binsz = 0.05

        # Keep track of the original value of edisp and switch-off energy
        # dispersion for the first iterations to speed-up the computations
        edisp_orig    = self['edisp'].boolean()
        self['edisp'] = False

        # Set flux ratio precision required for convergence to 5%
        ratio_precision = 0.05

        # Set energy boundaries
        self._set_obs_ebounds(emin, emax)

        # Determine mean energy for energy boundary
        e_mean   = math.sqrt(emin.TeV()*emax.TeV())
        loge     = math.log10(e_mean)
        erg_mean = e_mean * tev2erg

        # Compute Crab unit. This is the factor with which the Prefactor needs
        # to be multiplied to get 1 Crab.
        crab_flux = self._get_crab_flux(emin, emax)
        src_flux  = test_model[self._srcname].spectral().flux(emin, emax)
        crab_unit = crab_flux/src_flux

        # Initialise regression coefficient
        regcoeff = 0.0

        # Write header for energy bin
        self._log_string(gammalib.TERSE, '')
        self._log_header2(gammalib.TERSE, 'Energies: '+str(emin)+' - '+str(emax))

        # Write initial parameters
        self._log_header3(gammalib.TERSE, 'Initial parameters')
        self._log_value(gammalib.TERSE, 'Crab flux', str(crab_flux)+' ph/cm2/s')
        self._log_value(gammalib.TERSE, 'Source model flux', str(src_flux)+' ph/cm2/s')
        self._log_value(gammalib.TERSE, 'Crab unit factor', crab_unit)

        # Initialise loop
        results        = []
        iterations     = 0
        test_crab_flux = 0.1 # Initial test flux in Crab units (100 mCrab)

        # Write header for iterations for terse chatter level
        if self._logTerse():
            self._log_header3(gammalib.TERSE, 'Iterations')

        # Loop until we break
        while True:

            # Update iteration counter
            iterations += 1

            # Increment the seed value
            self._seed += 1

            # Write header for iteration into logger
            self._log_header2(gammalib.EXPLICIT, 'Iteration '+str(iterations))

            crab_prefactor = self._set_src_prefactor(test_model, crab_unit, test_crab_flux)

            # Simulate events for the models. "sim" holds an observation
            # container with observations containing the simulated events.
            sim = obsutils.sim(self.obs(), nbins=enumbins, seed=self._seed,
                               binsz=binsz, npix=npix,
                               log=self._log_clients,
                               debug=self['debug'].boolean(),
                               edisp=self['edisp'].boolean(),
                               nthreads=1)

            # Determine number of events in simulation
            nevents = sim.nobserved()

            # Write simulation results into logger
            self._log_header3(gammalib.EXPLICIT, 'Simulation')
            self._log_value(gammalib.EXPLICIT, 'Number of simulated events', nevents)

            # Fit test source to the simulated events in the observation
            # container
            fit = ctools.ctlike(sim)
            fit['edisp']    = self['edisp'].boolean()
            fit['nthreads'] = 1  # Avoids OpenMP conflict
            fit['debug']    = self['debug'].boolean()
            fit['chatter']  = self['chatter'].integer()
            fit.run()

            # Get model fitting results
            logL   = fit.opt().value()
            npred  = fit.obs().npred()
            models = fit.obs().models()
            source = models[self._srcname]
            ts     = source.ts()

            # Get fitted Crab, photon and energy fluxes
            prefactor   = modutils.normalisation_parameter(source)
            crab_flux   = prefactor.value() / crab_prefactor
            photon_flux = source.spectral().flux(emin, emax)
            energy_flux = source.spectral().eflux(emin, emax)

            # Compute differential sensitivity in unit erg/cm2/s by evaluating
            # the spectral model at the "e_mean" energy and by multipling the
            # result with the energy squared. Since the "eval()" method returns
            # an intensity in units of ph/cm2/s/MeV we multiply by 1.0e6 to
            # convert into ph/cm2/s/TeV, by "e_mean" to convert into ph/cm2/s,
            # and finally by "erg_mean" to convert to erg/cm2/s.
            energy      = gammalib.GEnergy(e_mean, 'TeV')
            sensitivity = source.spectral().eval(energy) * e_mean*erg_mean*1.0e6

            # Write fit results into logger
            name  = 'Iteration %d' % iterations
            value = ('TS=%10.4f  Sim=%9.4f mCrab  Fit=%9.4f mCrab  '
                     'Sens=%e erg/cm2/s' %
                     (ts, test_crab_flux*1000.0, crab_flux*1000.0, sensitivity))
            self._log_value(gammalib.TERSE, name, value)

            # If TS was non-positive then increase the test flux and start over
            if ts <= 0.0:

                # If the number of iterations was exceeded then stop
                if (iterations >= max_iter):
                    self._log_string(gammalib.TERSE,
                         ' Test ended after %d iterations.' % max_iter)
                    break

                # Increase test flux
                test_crab_flux *= 3.0

                # Signal start we start over
                self._log_string(gammalib.EXPLICIT,
                     'Non positive TS, increase test flux and start over.')

                # ... and start over
                continue

            # Append result entry to result list
            result = {'ts': ts, 'crab_flux': crab_flux,
                                'photon_flux': photon_flux,
                                'energy_flux': energy_flux}
            results.append(result)

            # Predict Crab flux at threshold TS using a linear regression of
            # the log(TS) and log(crab_flux) values that have so far been
            # computed. If not enough results are available then use a simple
            # TS scaling relation.
            correct = math.sqrt(ts_thres / ts)
            if len(results) > 1:
                try:
                    pred_crab_flux, regcoeff = self._predict_flux(results, ts_thres)
                    correct                  = pred_crab_flux / crab_flux
                except:
                    self._log_value(gammalib.TERSE, 'Skipping failed regression',
                                    'Retain simple TS scaling relation')

            # Compute extrapolated fluxes based on the flux correction factor
            crab_flux   = correct * crab_flux
            photon_flux = correct * photon_flux
            energy_flux = correct * energy_flux
            sensitivity = correct * sensitivity

            # If we have at least 3 results then check if the flux determination
            # at the TS threshold has converged
            if len(results) > 3:
                if test_crab_flux > 0:

                    # Compute fractional change in the Crab flux between two
                    # iterations
                    ratio = crab_flux/test_crab_flux

                    # If fractional change is smaller than the required position
                    # the iterations are stopped
                    if ratio > 1.0-ratio_precision and \
                       ratio < 1.0+ratio_precision:
                        value = ('TS=%10.4f  Sim=%9.4f mCrab                  '
                                 '     Sens=%e erg/cm2/s' %
                                 (ts, crab_flux*1000.0, sensitivity))
                        self._log_value(gammalib.TERSE, 'Converged result', value)
                        self._log_value(gammalib.TERSE, 'Converged flux ratio', ratio)
                        self._log_value(gammalib.TERSE, 'Regression coefficient',
                                        regcoeff)

                        # If the flux has converged then check if the original
                        # value of edisp was set, and is not what has fo far
                        # been used
                        if edisp_orig and not self['edisp'].boolean():

                            # Set edisp on and continue the calculation
                            self['edisp'] = True

                            # Log action
                            self._log_value(gammalib.TERSE, 'Converged result',
                                            'Now use energy dispersion after '
                                            'initial convergence without it')

                        # ... otherwise finish the calculation after convergence
                        else:
                            break

                # ... otherwise break with a zero flux
                else:
                    self._log_value(gammalib.TERSE, 'Not converged', 'Flux is zero')
                    break

            # Set test flux for next iteration
            test_crab_flux = crab_flux

            # Exit loop if number of trials exhausted
            if (iterations >= max_iter):
                self._log_string(gammalib.TERSE,
                                 ' Test ended after %d iterations.' % max_iter)
                break

        # Write header for iterations for terse chatter level
        if self._logTerse():
            self._log_header3(gammalib.TERSE, 'Iterations for source counts cuts')

        # Recover original energy dispersion setting
        self['edisp'] = edisp_orig

        # Take into account the excess event cuts
        n_bck_evts = None
        for iterations in range(max_iter):

            # Simulate event excess
            pass_evt_cut, n_bck_evts = self._sim_evt_excess(enumbins=enumbins,
                                                            binsz=binsz,
                                                            npix=npix,
                                                            n_bck_evts=n_bck_evts)

            # Write results into logger
            name  = 'Iteration %d' % iterations
            value = 'Fit=%9.4f mCrab  Sens=%e erg/cm2/s' % (crab_flux*1000.0, sensitivity)
            self._log_value(gammalib.TERSE, name, value)

            # If event excess cuts are passed then break now
            if pass_evt_cut:
                break

            # ... otherwise increase the test flux
            correct     = 1.0 + ratio_precision
            crab_flux   = correct * crab_flux
            photon_flux = correct * photon_flux
            energy_flux = correct * energy_flux
            sensitivity = correct * sensitivity

            # ...
            _ = self._set_src_prefactor(test_model, crab_unit, crab_flux)

        # Write fit results into logger
        self._log_header3(gammalib.TERSE, 'Fit results')
        self._log_value(gammalib.TERSE, 'Photon flux',
                        str(photon_flux)+' ph/cm2/s')
        self._log_value(gammalib.TERSE, 'Energy flux',
                        str(energy_flux)+' erg/cm2/s')
        self._log_value(gammalib.TERSE, 'Crab flux',
                        str(crab_flux*1000.0)+' mCrab')
        self._log_value(gammalib.TERSE, 'Differential sensitivity',
                        str(sensitivity)+' erg/cm2/s')
        self._log_value(gammalib.TERSE, 'Number of simulated events', nevents)
        self._log_header3(gammalib.TERSE, 'Test source model fitting')
        self._log_value(gammalib.TERSE, 'log likelihood', logL)
        self._log_value(gammalib.TERSE, 'Number of predicted events', npred)
        for model in models:
            self._log_value(gammalib.TERSE, 'Model', model.name())
            for par in model:
                self._log_string(gammalib.TERSE, str(par))

        # Restore energy boundaries of observation container
        for i, obs in enumerate(self.obs()):
            obs.events().ebounds(self._obs_ebounds[i])

        # Store result
        result = {'loge': loge, 'emin': emin.TeV(), 'emax': emax.TeV(), \
                  'crab_flux': crab_flux, 'photon_flux': photon_flux, \
                  'energy_flux': energy_flux, \
                  'sensitivity': sensitivity, 'regcoeff': regcoeff, \
                  'nevents': nevents, 'npred': npred}

        # Return result
        return result

    def _set_src_prefactor(self, test_model, crab_unit, test_crab_flux):
        """
        Create a copy of the test models, set the normalisation parameter
        of the test source in the models, and append the models to the observation.

        Parameters
        ----------
        test_model : `~gammalib.GModels`
            Test source model
        crab_unit : float
            Crab unit factor
        test_crab_flux : float
            Test flux in Crab units (100 mCrab)

        Returns
        -------
        crab_prefactor : float
            the prefactor that corresponds to a flux of 1 Crab.
        """
        # Initialise variables
        models         = test_model.copy()
        prefactor      = modutils.normalisation_parameter(models[self._srcname])
        crab_prefactor = prefactor.value() * crab_unit
        val_margin     = 0.01
        min_pref       = prefactor.min() * (1.0 + val_margin)
        max_pref       = prefactor.max() * (1.0 - val_margin)
        val_now        = crab_prefactor * test_crab_flux

        # Check whether prefactor is in valid range
        if val_now < min_pref or val_now > max_pref:

            # Store old prefactor value for logging
            val_old = val_now

            # Set new prefactor value
            val_now        = max(val_now, min_pref)
            val_now        = min(val_now, max_pref)
            crab_prefactor = val_now / test_crab_flux

            # Log prefactor modification
            self._log_value(gammalib.EXPLICIT, 'Prefactor range',
                            '['+str(min_pref)+','+str(max_pref)+']')
            self._log_value(gammalib.EXPLICIT, 'Initial Prefactor', val_old)
            self._log_value(gammalib.EXPLICIT, 'Updated Prefactor', val_now)

        # Store prefactor value
        prefactor.value(val_now)

        # Update the models
        self.obs().models(models)

        # Return Prefactor
        return crab_prefactor

    def _simulate_events(self, obs, rad, enumbins, binsz, npix):
        """
        Simulate events

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Observation container
        rad : float
            Simulation selection radius (degrees)
        enumbins : integer
            Number of energy bins
        binsz : float
            Pixel size
        npix : integer
            Number of pixels
        """
        # Initialise list of number of simulated events
        n_sim_evts = []

        # Simulate events
        for n_sim in range(15):

            # Increment the seed value
            self._seed += 1

            # Simulate observations
            sim = obsutils.sim(obs, nbins=enumbins, seed=self._seed, binsz=binsz,
                               npix=npix, log=self._log_clients,
                               debug=self['debug'].boolean(),
                               edisp=self['edisp'].boolean(), nthreads=1)

            # If a selection radius is given then select only the events within
            # that selection radius
            if rad > 0.0:

                # Run event selection tool
                select = ctools.ctselect(sim)
                select['rad']  = rad
                select['emin'] = 'INDEF'
                select['tmin'] = 'INDEF'
                select.run()

                # Store number of events
                n_sim_evts.append(select.obs().nobserved())

            # ... otherwise we simply store the number of simulated events
            else:
                n_sim_evts.append(sim.nobserved())

        # Determine median number of events
        median = self._median(n_sim_evts)

        # Log results
        if rad > 0.0:
            name = 'Median source events'
        else:
            name = 'Median background events'
        self._log_value(gammalib.EXPLICIT, name, median)

        # Return
        return median

    def _sim_evt_excess(self, enumbins, binsz, npix, n_bck_evts):
        """
        Return the number of excess events for the source model, compared
        to all other models

        Parameters
        ----------
        enumbins : int
            Number of bins for the observation simulation
        binsz : float
            Bin size for the observation simulation
        npix : int
            Pixel size for the observation simulation
        n_bck_evts : int
            Number of background counts from previous call

        Returns
        -------
        pass_evt_cut : bool
            Signals that the source passes the minimal excess criteria
        n_bck_evts : int
            Number of background counts
        """
        # Initialise results
        pass_evt_cut = True
        n_bck_evts   = 0

        # Get user parameters
        mincounts = self['mincounts'].integer()
        bkgexcess = self['bkgexcess'].real()
        bkgrad    = self['bkgrad'].real()

        # Continue only if cuts were specified
        if mincounts > 0 or bkgexcess > 0.0:

            # Get models and split them into source and background
            models     = self.obs().models()
            src_model  = gammalib.GModels()
            bck_models = gammalib.GModels()
            for model in models:
                if model.name() == self._srcname:
                    src_model.append(model)
                else:
                    bck_models.append(model)

            # Create observations for source and background
            src_obs = self.obs().copy()
            bck_obs = self.obs().copy()
            src_obs.models(src_model)
            bck_obs.models(bck_models)

            # Simulate source events
            n_src_evts = self._simulate_events(src_obs, 0.0, enumbins, binsz, npix)

            # Simulate background events in case that a background fraction cut
            # is applied and if the number of background events was not yet
            # estimated
            if bkgexcess > 0.0 and n_bck_evts == None:
                n_bck_evts = self._simulate_events(bck_obs, bkgrad, enumbins, binsz, npix)

            # Set cut results
            min_bkg_events = n_bck_evts * bkgexcess
            has_min_evts   = n_src_evts >= mincounts
            is_above_bck   = n_src_evts >= min_bkg_events
            pass_evt_cut   = has_min_evts and is_above_bck

            # Log results
            self._log_value(gammalib.EXPLICIT, 'Pass minimum counts cut', has_min_evts)
            self._log_value(gammalib.EXPLICIT, 'Pass background cut', is_above_bck)
            self._log_value(gammalib.EXPLICIT, 'Pass event cut', pass_evt_cut)
            self._log_value(gammalib.EXPLICIT, 'Minimum counts threshold', mincounts)
            self._log_value(gammalib.EXPLICIT, 'Background threshold', min_bkg_events)
            self._log_value(gammalib.EXPLICIT, 'Source events', n_src_evts)
            self._log_value(gammalib.EXPLICIT, 'Background events', n_bck_evts)

        # Return passing flag and number of background events
        return pass_evt_cut, n_bck_evts

    def _predict_flux(self, results, ts):
        """
        Predict Crab flux for a given TS value

        The Crab flux at a given Test Statistic value is predicted by doing a
        linear regression of the log(TS) and log(crab_flux) values in a results
        list.

        See https://en.wikipedia.org/wiki/Simple_linear_regression

        Parameters
        ----------
        results : list of dict
            List of results
        ts : float
            Test Statistic value

        Returns
        -------
        crab_flux_prediction, regcoeff : tuple of float
            Predicted Crab flux and regression coefficient
        """
        # Compute means and regression coefficient
        mean_x  = 0.0
        mean_y  = 0.0
        mean_xy = 0.0
        mean_xx = 0.0
        mean_yy = 0.0
        for result in results:
            x        = math.log(float(result['ts']))
            y        = math.log(float(result['crab_flux']))
            mean_x  += x
            mean_y  += y
            mean_xy += x * y
            mean_xx += x * x
            mean_yy += y * y
        norm     = 1.0 / float(len(results))
        mean_x  *= norm
        mean_y  *= norm
        mean_xy *= norm
        mean_xx *= norm
        mean_yy *= norm

        # Compute regression coefficient
        rxy_norm = (mean_xx - mean_x * mean_x) * (mean_yy - mean_y * mean_y)
        if rxy_norm < 1e-10:
            rxy_norm = 1
        else:
            rxy_norm = math.sqrt(rxy_norm)
        rxy      = (mean_xy - mean_x * mean_y) / rxy_norm
        regcoeff = rxy*rxy

        # Compute regression line slope
        beta_nom   = 0.0
        beta_denom = 0.0
        for result in results:
            x           = math.log(float(result['ts']))
            y           = math.log(float(result['crab_flux']))
            beta_nom   += (x - mean_x) * (y - mean_y)
            beta_denom += (x - mean_x) * (x - mean_x)
        beta = beta_nom / beta_denom

        # Compute regression line offset
        alpha = mean_y - beta * mean_x

        # Predict Crab flux at TS threshold
        log_ts_thres         = math.log(ts)
        crab_flux_prediction = math.exp(alpha + beta * log_ts_thres)

        # Return
        return (crab_flux_prediction, regcoeff)

    def _e_bin(self, ieng):
        """
        Determines sensivity in energy bin

        Parameters
        ----------
        ieng : int
            Energy bin number

        Returns
        -------
        result : dict
            Result dictionary
        """
        # Get sensitivity type
        sensitivity_type = self['type'].string()

        # Set energies
        if sensitivity_type == 'Differential':
            emin = self._ebounds.emin(ieng)
            emax = self._ebounds.emax(ieng)
        elif sensitivity_type == 'Integral':
            emin = self._ebounds.emin(ieng)
            emax = self._ebounds.emax()
        else:
            msg = ('Invalid sensitivity type "%s" encountered. Either '
                   'specify "Differential" or "Integral".' %
                   sensitivity_type)
            raise RuntimeError(msg)

        # Determine sensitivity
        result = self._get_sensitivity(emin, emax, self._models)

        # Return results
        return result

    def _create_fits(self, results):
        """
        Create FITS file from results

        Parameters
        ----------
        results : list of dict
            List of result dictionaries
        """
        # Create FITS table columns
        nrows        = len(results)
        e_mean       = gammalib.GFitsTableDoubleCol('E_MEAN', nrows)
        e_min        = gammalib.GFitsTableDoubleCol('E_MIN', nrows)
        e_max        = gammalib.GFitsTableDoubleCol('E_MAX', nrows)
        flux_crab    = gammalib.GFitsTableDoubleCol('FLUX_CRAB', nrows)
        flux_photon  = gammalib.GFitsTableDoubleCol('FLUX_PHOTON', nrows)
        flux_energy  = gammalib.GFitsTableDoubleCol('FLUX_ENERGY', nrows)
        sensitivity  = gammalib.GFitsTableDoubleCol('SENSITIVITY', nrows)
        regcoeff     = gammalib.GFitsTableDoubleCol('REGRESSION_COEFF', nrows)
        nevents      = gammalib.GFitsTableDoubleCol('NEVENTS', nrows)
        npred        = gammalib.GFitsTableDoubleCol('NPRED', nrows)
        e_mean.unit('TeV')
        e_min.unit('TeV')
        e_max.unit('TeV')
        flux_crab.unit('')
        flux_photon.unit('ph/cm2/s')
        flux_energy.unit('erg/cm2/s')
        sensitivity.unit('erg/cm2/s')
        regcoeff.unit('')
        nevents.unit('counts')
        npred.unit('counts')

        # Fill FITS table columns
        for i, result in enumerate(results):
            e_mean[i]      = math.pow(10.0, result['loge'])
            e_min[i]       = result['emin']
            e_max[i]       = result['emax']
            flux_crab[i]   = result['crab_flux']
            flux_photon[i] = result['photon_flux']
            flux_energy[i] = result['energy_flux']
            sensitivity[i] = result['sensitivity']
            regcoeff[i]    = result['regcoeff']
            nevents[i]     = result['nevents']
            npred[i]       = result['npred']

        # Initialise FITS Table with extension "SENSITIVITY"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('SENSITIVITY')

        # Add keywors for compatibility with gammalib.GMWLSpectrum
        table.card('INSTRUME', 'CTA', 'Name of Instrument')
        table.card('TELESCOP', 'CTA', 'Name of Telescope')

        # Stamp header
        self._stamp(table)

        # Add script keywords
        table.card('OBJECT',   self._srcname, 'Source for which sensitivity was estimated')
        table.card('TYPE',     self['type'].string(), 'Sensitivity type')
        table.card('SIGMA',    self['sigma'].real(), '[sigma] Sensitivity threshold')
        table.card('MAX_ITER', self['max_iter'].integer(), 'Maximum number of iterations')
        table.card('STAT',     self['statistic'].string(), 'Optimization statistic')
        table.card('MINCOUNT', self['mincounts'].integer(), 'Minimum number of source counts')
        table.card('BKGEXCES', self['bkgexcess'].real(), 'Background uncertainty fraction')
        table.card('BKGRAD',   self['bkgrad'].real(), '[deg] Background radius')
        table.card('SEED',     self['seed'].integer(), 'Seed value for random numbers')

        # Append filled columns to fits table
        table.append(e_mean)
        table.append(e_min)
        table.append(e_max)
        table.append(flux_crab)
        table.append(flux_photon)
        table.append(flux_energy)
        table.append(sensitivity)
        table.append(regcoeff)
        table.append(nevents)
        table.append(npred)

        # Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append(table)

        # Return
        return


    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Loop over observations and store a deep copy of the energy
        # boundaries for later use
        for obs in self.obs():
            self._obs_ebounds.append(obs.events().ebounds().copy())

        # Initialise results
        results  = []

        # Set test source model for this observation
        self._models = modutils.test_source(self.obs().models(), self._srcname,
                                            ra=self._ra, dec=self._dec)

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Input observation')

        # Write models into logger
        self._log_models(gammalib.NORMAL, self._models, 'Input model')

        # Write header
        self._log_header1(gammalib.TERSE, 'Sensitivity determination')
        self._log_value(gammalib.TERSE, 'Type', self['type'].string())

        # If using multiprocessing
        if self._nthreads > 1 and self._ebounds.size() > 1:

            # Compute energy bins
            args        = [(self, '_e_bin', i)
                           for i in range(self._ebounds.size())]
            poolresults = mputils.process(self._nthreads, mputils.mpfunc, args)

            # Construct results
            for ieng in range(self._ebounds.size()):
                results.append(poolresults[ieng][0])
                self._log_string(gammalib.TERSE, poolresults[ieng][1]['log'], False)

        # Otherwise, loop over energy bins
        else:
            for ieng in range(self._ebounds.size()):

                # Run analysis in energy bin
                result = self._e_bin(ieng)

                # Append results
                results.append(result)

        # Create FITS file
        self._create_fits(results)

        # Return
        return

    def save(self):
        """
        Save sensitivity FITS file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save sensitivity curve')

        # Continue only if FITS file is valid
        if self._fits != None:

            # Get outmap parameter
            outfile = self['outfile'].filename()

            # Log file name
            self._log_value(gammalib.NORMAL, 'Sensitivity file', outfile.url())

            # Save sensitivity
            self._fits.saveto(outfile, self['clobber'].boolean())

        # Return
        return

    def sensitivity(self):
        """
        Return sensitivity FITS file

        Returns:
            FITS file containing sensitivity curve
        """
        # Return
        return self._fits


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cssens(sys.argv)

    # Run application
    app.execute()
