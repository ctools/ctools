#! /usr/bin/env python
# ==========================================================================
# This script computes the CTA sensitivity.
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
        self._debug       = False
        
        # Initialise application by calling the appropriate class
        # constructor.
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

        # Set some fixed parameters
        self._debug = self['debug'].boolean() # Debugging in client tools

        # Derive some parameters
        self._ebounds = gammalib.GEbounds(self._bins,
                                          gammalib.GEnergy(self._emin, 'TeV'),
                                          gammalib.GEnergy(self._emax, 'TeV'))

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log('\n')

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

    def _set_models(self, fitspat=False, fitspec=False):
        """
        Set full model and background model

        Parameters
        ----------
        fitspat : bool, optional
            Fit spatial parameter?
        fitspec : bool, optional
            Fit spectral parameters?

        Returns
        -------
        full_model, bkg_model : `~gammalib.GModels`
            Full model and background model containers
        """
        # Retrieve full model from observation container
        full_model = self._obs.models().copy()

        # Get source model
        model = full_model[self._srcname]

        # If source model has no "Prefactor" parameter then raise an exception
        if not model.has_par('Prefactor'):
            msg = ('Model "%s" has no parameter "Prefactor". Only spectral '
                   'models with a "Prefactor" parameter are supported.' %
                   self._srcname)
            raise RuntimeError(msg)

        # Set source position
        if self._ra != None and self._dec != None:
            if model.has_par('RA') and model.has_par('DEC'):
                model['RA'].value(self._ra)
                model['DEC'].value(self._dec)

        # Set possible spatial and spectral parameters
        spatial  = ['RA', 'DEC', 'Sigma', 'Radius', 'Width', 'PA',
                    'MinorRadius', 'MajorRadius']
        spectral = ['Index', 'Index1', 'Index2', 'BreakEnergy', 'CutoffEnergy',
                    'InverseCutoffEnergy']

        # Fit or fix spatial parameters
        for par in spatial:
            if model.has_par(par):
                if fitspat:
                    model[par].free()
                else:
                    model[par].fix()

        # Fit or fix spectral parameters
        for par in spectral:
            if model.has_par(par):
                if fitspec:
                    model[par].free()
                else:
                    model[par].fix()

        # Create background model
        bkg_model = full_model.copy()
        bkg_model.remove(self._srcname)

        # Return models
        return full_model, bkg_model

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
        for obs in self._obs:

            # Get observation energy boundaries
            obs_ebounds = obs.events().ebounds()
            
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

    def _get_sensitivity(self, emin, emax, bkg_model, full_model):
        """
        Determine sensitivity for given observations.

        Parameters
        ----------
        emin : `~gammalib.GEnergy`
            Minimum energy for fitting and flux computation
        emax : `~gammalib.GEnergy`
            Maximum energy for fitting and flux computation
        bkg_model : `~gammalib.GModels`
            Background model
        full_model : `~gammalib.GModels`
            Source plus background model

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

        # Compute Crab unit (this is the factor with which the Prefactor needs
        # to be multiplied to get 1 Crab
        crab_flux = self._get_crab_flux(emin, emax)
        src_flux  = full_model[self._srcname].spectral().flux(emin, emax)
        crab_unit = crab_flux/src_flux

        # Write header
        if self._logTerse():
            self._log('\n')
            self._log.header2('Energies: '+str(emin)+' - '+str(emax))
            self._log.parformat('Crab flux')
            self._log(crab_flux)
            self._log(' ph/cm2/s\n')
            self._log.parformat('Source model flux')
            self._log(src_flux)
            self._log(' ph/cm2/s\n')
            self._log.parformat('Crab unit factor')
            self._log(crab_unit)
            self._log('\n')

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

            # Write header
            if self._logExplicit():
                self._log.header2('Iteration '+str(iterations))

            # Set source model. crab_prefactor is the Prefactor that
            # corresponds to 1 Crab
            src_model      = full_model.copy()
            crab_prefactor = src_model[self._srcname]['Prefactor'].value() * \
                             crab_unit
            src_model[self._srcname]['Prefactor'].value(crab_prefactor * \
                             test_crab_flux)
            self._obs.models(src_model)

            # Simulate events
            sim = obsutils.sim(self._obs, nbins=self._enumbins, seed=iterations,
                               binsz=self._binsz, npix=self._npix,
                               log=self._log_clients, debug=self._debug,
                               edisp=self._edisp)

            # Determine number of events in simulation
            nevents = 0.0
            for run in sim:
                nevents += run.events().number()

            # Write simulation results
            if self._logExplicit():
                self._log.header3('Simulation')
                self._log.parformat('Number of simulated events')
                self._log(nevents)
                self._log('\n')

            # Fit background only
            sim.models(bkg_model)
            like       = obsutils.fit(sim, log=self._log_clients,
                                      debug=self._debug, edisp=self._edisp)
            result_bgm = like.obs().models().copy()
            LogL_bgm   = like.opt().value()
            npred_bgm  = like.obs().npred()

            # Assess quality based on a comparison between Npred and Nevents
            quality_bgm = npred_bgm-nevents

            # Write background fit results
            if self._logExplicit():
                self._log.header3('Background model fit')
                self._log.parformat('log likelihood')
                self._log(LogL_bgm)
                self._log('\n')
                self._log.parformat('Number of predicted events')
                self._log(npred_bgm)
                self._log('\n')
                self._log.parformat('Fit quality')
                self._log(quality_bgm)
                self._log('\n')

            # Write model fit results
            if self._logExplicit():
                for model in result_bgm:
                    self._log.parformat('Model')
                    self._log(model.name())
                    self._log('\n')
                    for par in model:
                        self._log(str(par)+'\n')

            # Fit background and test source
            sim.models(src_model)
            like       = obsutils.fit(sim, log=self._log_clients,
                                      debug=self._debug, edisp=self._edisp)
            result_all = like.obs().models().copy()
            LogL_all   = like.opt().value()
            npred_all  = like.obs().npred()
            ts         = 2.0*(LogL_bgm-LogL_all)

            # Assess quality based on a comparison between Npred and Nevents
            quality_all = npred_all-nevents

            # Write background and test source fit results
            if self._logExplicit():
                self._log.header3('Background and test source model fit')
                self._log.parformat('Test statistics')
                self._log(ts)
                self._log('\n')
                self._log.parformat('log likelihood')
                self._log(LogL_all)
                self._log('\n')
                self._log.parformat('Number of predicted events')
                self._log(npred_all)
                self._log('\n')
                self._log.parformat('Fit quality')
                self._log(quality_all)
                self._log('\n')
                #
                for model in result_all:
                    self._log.parformat('Model')
                    self._log(model.name())
                    self._log('\n')
                    for par in model:
                        self._log(str(par)+'\n')

            # Start over if TS was non-positive
            if ts <= 0.0:
                if self._logExplicit():
                    self._log('Non positive TS. Start over.\n')
                continue

            # Get fitted Crab, photon and energy fluxes
            crab_flux     = result_all[self._srcname]['Prefactor'].value() / \
                            crab_prefactor
            #crab_flux_err = result_all[self._srcname]['Prefactor'].error() / \
            #                crab_prefactor
            photon_flux   = result_all[self._srcname].spectral().flux(emin, emax)
            energy_flux   = result_all[self._srcname].spectral().eflux(emin, emax)

            # Compute differential sensitivity in unit erg/cm2/s
            energy      = gammalib.GEnergy(e_mean, 'TeV')
            time        = gammalib.GTime()
            sensitivity = result_all[self._srcname].spectral().eval(energy, time) * \
                          e_mean*erg_mean * 1.0e6

            # Compute flux correction factor based on average TS
            correct = 1.0
            if ts > 0:
                correct = math.sqrt(self._ts_thres/ts)

            # Compute extrapolated fluxes
            crab_flux   = correct * crab_flux
            photon_flux = correct * photon_flux
            energy_flux = correct * energy_flux
            sensitivity = correct * sensitivity
            crab_flux_value.append(crab_flux)
            photon_flux_value.append(photon_flux)
            energy_flux_value.append(energy_flux)
            sensitivity_value.append(sensitivity)

            # Write background and test source fit results
            if self._logExplicit():
                self._log.parformat('Photon flux')
                self._log(photon_flux)
                self._log(' ph/cm2/s\n')
                self._log.parformat('Energy flux')
                self._log(energy_flux)
                self._log(' erg/cm2/s\n')
                self._log.parformat('Crab flux')
                self._log(crab_flux*1000.0)
                self._log(' mCrab\n')
                self._log.parformat('Differential sensitivity')
                self._log(sensitivity)
                self._log(' erg/cm2/s\n')
                for model in result_all:
                    self._log.parformat('Model')
                    self._log(model.name())
                    self._log('\n')
                    for par in model:
                        self._log(str(par)+'\n')
            elif self._logTerse():
                self._log.parformat('Iteration '+str(iterations))
                self._log('TS=')
                self._log(ts)
                self._log(' ')
                self._log('corr=')
                self._log(correct)
                self._log('  ')
                self._log(photon_flux)
                self._log(' ph/cm2/s = ')
                self._log(energy_flux)
                self._log(' erg/cm2/s = ')
                self._log(crab_flux*1000.0)
                self._log(' mCrab = ')
                self._log(sensitivity)
                self._log(' erg/cm2/s\n')

            # Compute sliding average of extrapolated fitted prefactor,
            # photon and energy flux. This damps out fluctuations and
            # improves convergence
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

        # Write fit results
        if self._logTerse():
            self._log.header3('Fit results')
            self._log.parformat('Test statistics')
            self._log(ts)
            self._log('\n')
            self._log.parformat('Photon flux')
            self._log(photon_flux)
            self._log(' ph/cm2/s\n')
            self._log.parformat('Energy flux')
            self._log(energy_flux)
            self._log(' erg/cm2/s\n')
            self._log.parformat('Crab flux')
            self._log(crab_flux*1000.0)
            self._log(' mCrab\n')
            self._log.parformat('Differential sensitivity')
            self._log(sensitivity)
            self._log(' erg/cm2/s\n')
            self._log.parformat('Number of simulated events')
            self._log(nevents)
            self._log('\n')
            self._log.header3('Background and test source model fitting')
            self._log.parformat('log likelihood')
            self._log(LogL_all)
            self._log('\n')
            self._log.parformat('Number of predicted events')
            self._log(npred_all)
            self._log('\n')
            for model in result_all:
                self._log.parformat('Model')
                self._log(model.name())
                self._log('\n')
                for par in model:
                    self._log(str(par)+'\n')
            self._log.header3('Background model fit')
            self._log.parformat('log likelihood')
            self._log(LogL_bgm)
            self._log('\n')
            self._log.parformat('Number of predicted events')
            self._log(npred_bgm)
            self._log('\n')
            for model in result_bgm:
                self._log.parformat('Model')
                self._log(model.name())
                self._log('\n')
                for par in model:
                    self._log(str(par)+'\n')

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
        
        # Loop over observations and store ebounds
        for obs in self._obs:
            self._obs_ebounds.append(obs.events().ebounds())
        
        # Initialise script
        colnames = ['loge', 'emin', 'emax', 'crab_flux', 'photon_flux',
                    'energy_flux', 'sensitivity']
        results  = []

        # Initialise models
        full_model, bkg_model = self._set_models()

        # Write models into logger
        if self._logTerse():
            self._log('\n')
            self._log.header1('Models')
            self._log.header2('Background model')
            self._log(str(bkg_model))
            self._log('\n\n')
            self._log.header2('Full model')
            self._log(str(full_model))
            self._log('\n')

        # Write header
        if self._logTerse():
            self._log('\n')
            self._log.header1('Sensitivity determination')
            self._log.parformat('Type')
            self._log(self._type)
            self._log('\n')

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
            result = self._get_sensitivity(emin, emax, bkg_model, full_model)

            # Write results
            if ieng == 0:
                f      = open(self._outfile.url(), 'w')
                writer = csv.DictWriter(f, colnames)
                headers = {}
                for n in colnames:
                    headers[n] = n
                writer.writerow(headers)
                writer.writerow(result)
                f.close()
            else:
                f      = open(self._outfile.url(), 'a')
                writer = csv.DictWriter(f, colnames)
                writer.writerow(result)
                f.close()

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
