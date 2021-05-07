#! /usr/bin/env python
# ==========================================================================
# Computes the array sensitivity using the Test Statistic for a test source
#
# Copyright (C) 2011-2021 Juergen Knoedlseder
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
from cscripts import ioutils
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
        Constructor.
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise class members
        self._ebounds     = gammalib.GEbounds()
        self._obs_ebounds = []
        self._srcname     = ''
        self._ra          = None
        self._dec         = None
        self._log_clients = False
        self._models      = gammalib.GModels()
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
                 'srcname'      : self._srcname,
                 'ra'           : self._ra,
                 'dec'          : self._dec,
                 'log_clients'  : self._log_clients,
                 'models'       : self._models,
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
        self._srcname     = state['srcname']
        self._ra          = state['ra']
        self._dec         = state['dec']
        self._log_clients = state['log_clients']
        self._models      = state['models']
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

        # Get source name
        self._srcname = self['srcname'].string()

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
        self['outfile'].filename()
        self['edisp'].boolean()
        self['debug'].boolean()

        # Derive some parameters
        self._ebounds = gammalib.GEbounds(bins,
                                          gammalib.GEnergy(emin, 'TeV'),
                                          gammalib.GEnergy(emax, 'TeV'))

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)


        # Set number of processes for multiprocessing
        self._nthreads = mputils.nthreads(self)

        # Return
        return

    def _set_obs(self, emin, emax, lpnt=None, bpnt=None):
        """
        Set an observation container

        Parameters
        ----------
        emin : float
            Minimum energy (TeV)
        emax : float
            Maximum energy (TeV)
        lpnt : float, optional
            Galactic longitude of pointing (deg)
        bpnt : float, optional
            Galactic latitude of pointing (deg)

        Returns
        -------
        obs : `~gammalib.GObservations`
            Observation container
        """
        # If an observation was provided on input then load it from XML file
        if self['inobs'].is_valid():
            obs = self._get_observations()

        # ... otherwise allocate a single observation
        else:

            # If no pointing is specified and if the test source has a position
            # then set the pointing to the position of the test source
            if lpnt is None and bpnt is None:
                models = gammalib.GModels(self['inmodel'].filename())
                source = models[self['srcname'].string()]
                if source.has_par('RA') and source.has_par('DEC'):
                    pntdir = gammalib.GSkyDir()
                    pntdir.radec_deg(source['RA'].value(), source['DEC'].value())
                    lpnt   = pntdir.l_deg()
                    bpnt   = pntdir.b_deg()
                elif source.has_par('GLON') and source.has_par('GLAT'):
                    lpnt = source['GLON'].value()
                    bpnt = source['GLAT'].value()
                else:
                    lpnt = 0.0
                    bpnt = 0.0

            # Set source position
            srcdir    = gammalib.GSkyDir()
            srcdir.lb_deg(lpnt, bpnt)
            self._ra  = srcdir.ra_deg()
            self._dec = srcdir.dec_deg()

            # Set pointing direction offset in galactic latitude
            offset  = self['offset'].real()
            bpnt   += offset
            if bpnt > 90.0:
                bpnt  = 180.0 - bpnt
                lpnt += 180.0
            pntdir  = gammalib.GSkyDir()
            pntdir.lb_deg(lpnt, bpnt)

            # Read relevant user parameters
            caldb    = self['caldb'].string()
            irf      = self['irf'].string()
            deadc    = self['deadc'].real()
            duration = self['duration'].real()
            rad      = self['rad'].real()

            # Allocate observation container
            obs = gammalib.GObservations()

            # Create CTA observation
            run = obsutils.set_obs(pntdir, caldb=caldb, irf=irf,
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

            # Write header for iteration into logger
            self._log_header2(gammalib.EXPLICIT, 'Iteration '+str(iterations))

            # Create a copy of the test models, set the normalisation parameter
            # of the test source in the models, and append the models to the
            # observation. "crab_prefactor" is the prefactor that corresponds
            # to a flux of 1 Crab.
            models         = test_model.copy()
            prefactor      = modutils.normalisation_parameter(models[self._srcname])
            crab_prefactor = prefactor.value() * crab_unit
            prefactor.value(crab_prefactor * test_crab_flux)
            self.obs().models(models)

            # Simulate events for the models. "sim" holds an observation
            # container with observations containing the simulated events.
            sim = obsutils.sim(self.obs(), nbins=enumbins, seed=iterations,
                               binsz=binsz, npix=npix,
                               log=self._log_clients,
                               debug=self['debug'].boolean(),
                               edisp=self['edisp'].boolean(),
                               nthreads=1)

            # Determine number of events in simulation by summing the events
            # over all observations in the observation container
            nevents = 0.0
            for run in sim:
                nevents += run.events().number()

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
            # computed. If not enough results are available than use a simple
            # TS scaling relation.
            if len(results) > 1:
                pred_crab_flux, regcoeff = self._predict_flux(results, ts_thres)
                correct                  = pred_crab_flux / crab_flux
            else:
                correct = math.sqrt(ts_thres/ts)

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
                        break
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
        rxy      = (mean_xy - mean_x * mean_y) / \
                   math.sqrt((mean_xx - mean_x*mean_x) *
                             (mean_yy - mean_y*mean_y))
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


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Loop over observations and store a deep copy of the energy
        # boundaries for later use
        for obs in self.obs():
            self._obs_ebounds.append(obs.events().ebounds().copy())

        # Initialise script
        colnames = ['loge', 'emin', 'emax', 'crab_flux', 'photon_flux',
                    'energy_flux', 'sensitivity', 'regcoeff', 'nevents',
                    'npred']
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

                #Run analysis in energy bin
                result = self._e_bin(ieng)

                # Append results
                results.append(result)

        # Write out trial result
        for ieng, result in enumerate(results):
            ioutils.write_csv_row(self['outfile'].filename().url(), ieng,
                                  colnames, result)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cssens(sys.argv)

    # Run application
    app.execute()
