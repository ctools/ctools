#! /usr/bin/env python
# ==========================================================================
# Computes the array sensitivity using a test source
#
# Copyright (C) 2011-2016 Juergen Knoedlseder
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
import csv
import math
import gammalib
import ctools
from cscripts import obsutils
from cscripts import modutils
from cscripts import ioutils


# ============ #
# cssens class #
# ============ #
class cssens(ctools.cscript):
    """
    Computes the CTA sensitivity

    This class computes the CTA sensitivity for a number of energy bins using
    ctlike. Spectra are fitted in narrow energy bins to simulated data,
    and the flux level is determined that leads to a particular significance.

    The significance is determined using the Test Statistic, defined as twice
    the likelihood difference between fitting with and  without the test source.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = 'cssens'
        self._version = '1.2.0'

        # Initialise class members
        self._obs         = gammalib.GObservations()
        self._ebounds     = gammalib.GEbounds()
        self._obs_ebounds = []
        self._srcname     = ''
        self._outfile     = gammalib.GFilename()
        self._ra          = None
        self._dec         = None
        self._edisp       = False
        self._emin        = 0.020
        self._emax        = 200.0
        self._bins        = 21
        self._enumbins    = 0
        self._npix        = 200
        self._binsz       = 0.05
        self._type        = 'Differential'
        self._ts_thres    = 25.0
        self._max_iter    = 50
        self._num_avg     = 3
        self._log_clients = False
        
        # Initialise application by calling the appropriate base class
        # constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get user parameters from parfile
        """
        # Set observation if not done before
        if self._obs == None or self._obs.size() == 0:
            self._obs = self._set_obs()

        # Set models if we have none
        if self._obs.models().size() == 0:
            self._obs.models(self['inmodel'].filename())

        # Get source name
        self._srcname = self['srcname'].string()

        # Read further parameters
        self._outfile = self['outfile'].filename()
        self._emin    = self['emin'].real()
        self._emax    = self['emax'].real()
        self._bins    = self['bins'].integer()

        # Read parameters for binned if requested
        self._enumbins = self['enumbins'].integer()
        if not self._enumbins == 0:
            self._npix  = self['npix'].integer()
            self._binsz = self['binsz'].real()

        # Read remaining parameters
        self._edisp    = self['edisp'].boolean()
        self._ts_thres = self['sigma'].real() * self['sigma'].real()
        self._max_iter = self['max_iter'].integer()
        self._num_avg  = self['num_avg'].integer()
        self._type     = self['type'].string()

        # Query remaining parameters
        self['debug'].boolean()

        # Derive some parameters
        self._ebounds = gammalib.GEbounds(self._bins,
                                          gammalib.GEnergy(self._emin, 'TeV'),
                                          gammalib.GEnergy(self._emax, 'TeV'))

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _set_obs(self, lpnt=0.0, bpnt=0.0, emin=0.1, emax=100.0):
        """
        Set an observation container

        Parameters
        ----------
        lpnt : float, optional
            Galactic longitude of pointing (deg)
        bpnt : float, optional
            Galactic latitude of pointing (deg)
        emin : float, optional
            Minimum energy (TeV)
        emax : float, optional
            Maximum energy (TeV)

        Returns
        -------
        obs : `~gammalib.GObservations`
            Observation container
        """
        # If an observation was provided on input then load it from XML file
        filename = self['inobs'].filename()
        if filename != 'NONE' and filename != '':
            obs = self._get_observations()

        # ... otherwise allocate a single observation
        else:

            # Read relevant user parameters
            caldb    = self['caldb'].string()
            irf      = self['irf'].string()
            deadc    = self['deadc'].real()
            duration = self['duration'].real()
            rad      = self['rad'].real()

            # Allocate observation container
            obs = gammalib.GObservations()

            # Set single pointing
            pntdir = gammalib.GSkyDir()
            pntdir.lb_deg(lpnt, bpnt)

            # Create CTA observation
            run = obsutils.set_obs(pntdir, caldb=caldb, irf=irf,
                                   duration=duration, deadc=deadc,
                                   emin=emin, emax=emax, rad=rad)

            # Append observation to container
            obs.append(run)

            # Set source position
            offset    = self['offset'].real()
            pntdir.lb_deg(lpnt, bpnt+offset)
            self._ra  = pntdir.ra_deg()
            self._dec = pntdir.dec_deg()

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
        for i, obs in enumerate(self._obs):

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

        # Set energy boundaries
        self._set_obs_ebounds(emin, emax)

        # Determine energy boundaries from first observation in the container
        loge      = math.log10(math.sqrt(emin.TeV()*emax.TeV()))
        e_mean    = math.pow(10.0, loge)
        erg_mean  = e_mean * tev2erg

        # Compute Crab unit. This is the factor with which the Prefactor needs
        # to be multiplied to get 1 Crab.
        crab_flux = self._get_crab_flux(emin, emax)
        src_flux  = test_model[self._srcname].spectral().flux(emin, emax)
        crab_unit = crab_flux/src_flux

        # Write header and some introductory information into logger
        self._log_header2(gammalib.TERSE, 'Energies: '+str(emin)+' - '+str(emax))
        self._log_value(gammalib.TERSE, 'Crab flux', str(crab_flux)+' ph/cm2/s')
        self._log_value(gammalib.TERSE, 'Source model flux', str(src_flux)+
                        ' ph/cm2/s')
        self._log_value(gammalib.TERSE, 'Crab unit factor', crab_unit)

        # Initialise loop
        crab_flux_value   = []
        photon_flux_value = []
        energy_flux_value = []
        sensitivity_value = []
        iterations        = 0
        test_crab_flux    = 0.1 # Initial test flux in Crab units (100 mCrab)

        # Loop until we break
        while True:

            # Update iteration counter
            iterations += 1

            # Write header for iteration into logger
            self._log_header2(gammalib.EXPLICIT, 'Iteration '+str(iterations))

            # Create a copy of the test models, set the prefactor of the test
            # source in the models, and append the models to the observation.
            # "crab_prefactor" is the Prefactor that corresponds to a flux of
            # 1 Crab.
            models         = test_model.copy()
            crab_prefactor = models[self._srcname]['Prefactor'].value() * crab_unit
            models[self._srcname]['Prefactor'].value(crab_prefactor * test_crab_flux)
            self._obs.models(models)

            # Simulate events for the models. "sim" holds an observation
            # container with observations containing the simulated events.
            sim = obsutils.sim(self._obs, nbins=self._enumbins, seed=iterations,
                               binsz=self._binsz, npix=self._npix,
                               log=self._log_clients,
                               debug=self['debug'].boolean(),
                               edisp=self._edisp)

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
            fit['edisp']   = edisp=self._edisp
            fit['debug']   = self['debug'].boolean()
            fit['chatter'] = self['chatter'].integer()
            fit.run()

            # Get model fitting results
            logL   = fit.opt().value()
            npred  = fit.obs().npred()
            models = fit.obs().models()
            source = models[self._srcname]
            ts     = source.ts()

            # Get fitted Crab, photon and energy fluxes
            crab_flux   = source['Prefactor'].value() / crab_prefactor
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

            # Assess quality based on a comparison between Npred and Nevents
            quality = npred - nevents

            # Write test source fit results into logger
            self._log_header3(gammalib.EXPLICIT, 'Test source model fit')
            self._log_value(gammalib.EXPLICIT, 'Test statistics', ts)
            self._log_value(gammalib.EXPLICIT, 'log likelihood', logL)
            self._log_value(gammalib.EXPLICIT, 'Number of predicted events', quality)
            self._log_value(gammalib.EXPLICIT, 'Fit quality', npred)
            for model in models:
                self._log_value(gammalib.EXPLICIT, 'Model', model.name())
                for par in model:
                    self._log_string(gammalib.EXPLICIT, str(par))

            # If TS was non-positive then increase the test flux and start
            # over
            if ts <= 0.0:

                # If the number of iterations was exceeded then stop
                if (iterations >= self._max_iter):
                    self._log_string(gammalib.TERSE,
                         ' Test ended after %d iterations.' % self._max_iter)
                    break

                # Increase test flux by a factor of 2
                test_crab_flux = test_crab_flux * 2.0

                # Signal start we start over
                self._log_string(gammalib.EXPLICIT,
                     'Non positive TS, increase test flux and start over.')

                # ... and start over
                continue

            # Compute flux correction factor based on average TS. We
            # assume here that the significance is proportional to the source
            # flux, and that the significance scales as the square root of
            # the Test Statistic
            correct = math.sqrt(self._ts_thres/ts)

            # Compute extrapolated fluxes based on the flux correction factor
            crab_flux   = correct * crab_flux
            photon_flux = correct * photon_flux
            energy_flux = correct * energy_flux
            sensitivity = correct * sensitivity
            crab_flux_value.append(crab_flux)
            photon_flux_value.append(photon_flux)
            energy_flux_value.append(energy_flux)
            sensitivity_value.append(sensitivity)

            # Write fit results
            if self._logExplicit():
                self._log_value(gammalib.EXPLICIT, 'Photon flux',
                                str(photon_flux)+' ph/cm2/s')
                self._log_value(gammalib.EXPLICIT, 'Energy flux',
                                str(energy_flux)+' erg/cm2/s')
                self._log_value(gammalib.EXPLICIT, 'Crab flux',
                                str(crab_flux*1000.0)+' mCrab')
                self._log_value(gammalib.EXPLICIT, 'Differential sensitivity',
                                str(sensitivity)+' erg/cm2/s')
                for model in models:
                    self._log_value(gammalib.EXPLICIT, 'Model', model.name())
                    for par in model:
                        self._log_string(gammalib.EXPLICIT, str(par))
            elif self._logTerse():
                name  = 'Iteration %d' % iterations
                value = ('TS=%.4f corr=%e  %e ph/cm2/s = %e erg/cm2/s = %.2f '
                         'mCrab = %e erg/cm2/s' % (ts, correct, photon_flux,
                         energy_flux, crab_flux*1000.0, sensitivity))
                self._log_value(gammalib.TERSE, name, value)

            # Compute sliding average of extrapolated fitted prefactor,
            # photon and energy flux. This damps out fluctuations and
            # improves convergence.
            crab_flux   = 0.0
            photon_flux = 0.0
            energy_flux = 0.0
            sensitivity = 0.0
            num         = 0.0
            for k in range(self._num_avg):
                inx = len(crab_flux_value) - k - 1
                if inx >= 0:
                    crab_flux   += crab_flux_value[inx]
                    photon_flux += photon_flux_value[inx]
                    energy_flux += energy_flux_value[inx]
                    sensitivity += sensitivity_value[inx]
                    num      += 1.0
            crab_flux   /= num
            photon_flux /= num
            energy_flux /= num
            sensitivity /= num

            # Compare average flux to last average
            if iterations > self._num_avg:
                if test_crab_flux > 0:
                    ratio = crab_flux/test_crab_flux

                    # We have 2 convergence criteria:
                    # 1. The average flux does not change
                    # 2. The flux correction factor is small
                    if ratio   >= 0.99 and ratio   <= 1.01 and \
                       correct >= 0.9  and correct <= 1.1:
                        if self._logTerse():
                            self._log(' Converged ('+str(ratio)+')\n')
                        break
                else:
                    if self._logTerse():
                        self._log(' Flux is zero.\n')
                    break

            # Use average for next iteration
            test_crab_flux = crab_flux

            # Exit loop if number of trials exhausted
            if (iterations >= self._max_iter):
                if self._logTerse():
                    self._log(' Test ended after '+str(self._max_iter)+
                              ' iterations.\n')
                break

        # Write fit results into logger
        self._log_header3(gammalib.TERSE, 'Fit results')
        self._log_value(gammalib.TERSE, 'Test statistics', ts)
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
        for i, obs in enumerate(self._obs):
            obs.events().ebounds(self._obs_ebounds[i])
        
        # Store result
        result = {'loge': loge, 'emin': emin.TeV(), 'emax': emax.TeV(), \
                  'crab_flux': crab_flux, 'photon_flux': photon_flux, \
                  'energy_flux': energy_flux, \
                  'sensitivity': sensitivity}

        # Return result
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
        for obs in self._obs:
            self._obs_ebounds.append(obs.events().ebounds().copy())
        
        # Initialise script
        colnames = ['loge', 'emin', 'emax', 'crab_flux', 'photon_flux',
                    'energy_flux', 'sensitivity']
        results  = []

        # Set test source model for this observation
        models = modutils.test_source(self._obs.models(), self._srcname,
                                      ra=self._ra, dec=self._dec)

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self._obs, 'Observation')

        # Write models into logger
        self._log_models(gammalib.NORMAL, models, 'Model')

        # Write header
        self._log_header1(gammalib.TERSE, 'Sensitivity determination')
        self._log_value(gammalib.TERSE, 'Type', str(self._type))

        # Loop over energy bins
        for ieng in range(self._ebounds.size()):

            # Set energies
            if self._type == 'Differential':
                emin  = self._ebounds.emin(ieng)
                emax  = self._ebounds.emax(ieng)
            elif self._type == 'Integral':
                emin  = self._ebounds.emin(ieng)
                emax  = self._ebounds.emax()
            else:
                msg = ('Invalid sensitivity type "%s" encountered. Either '
                       'specify "Differential" or "Integral".' % self._type)
                raise RuntimeError(msg)

            # Determine sensitivity
            result = self._get_sensitivity(emin, emax, models)

            # Write out trial result
            ioutils.write_csv_row(self._outfile.url(), ieng, colnames, result)

            # Append results
            results.append(result)

        # Return
        return

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Run the script
        self.run()

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
