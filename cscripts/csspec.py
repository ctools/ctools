#! /usr/bin/env python
# ==========================================================================
# Generates a spectrum.
#
# Copyright (C) 2014-2017 Michael Mayer
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
import gammalib
import ctools


# ============ #
# csspec class #
# ============ #
class csspec(ctools.csobservation):
    """
    Generates a spectrum

    This class implements the generation of a Spectral Energy Distribution
    (SED) from gamma-ray observations.
    """

    # Constructors and destructors
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise data members
        self._ebounds     = gammalib.GEbounds()
        self._fits        = None
        self._binned_mode = False
        self._onoff_mode  = False
        self._method      = 'AUTO'

        # Return
        return

    def __del__(self):
        """
        Destructor
        """
        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Set observation if not done before
        if self.obs().is_empty():
            self._require_inobs('csspec::get_parameters()')
            self.obs(self._get_observations())

        # Set observation statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

        # Set models if we have none
        if self.obs().models().is_empty():
            self.obs().models(self['inmodel'].filename())

        # Query source name
        self['srcname'].string()

        # Get spectrum generation method
        self._method = self['method'].string()

        # Collect number of unbinned, binned and On/Off observations in
        # observation container
        n_unbinned = 0
        n_binned   = 0
        n_onoff    = 0
        for obs in self.obs():
            if obs.classname() == 'GCTAObservation':
                if obs.eventtype() == 'CountsCube':
                    n_binned += 1
                else:
                    n_unbinned += 1
            elif obs.classname() == 'GCTAOnOffObservation':
                n_onoff += 1
        n_cta   = n_unbinned + n_binned + n_onoff
        n_other = self.obs().size() - n_cta

        # If spectrum method is not "NODES" then set spectrum method and
        # script mode according to type of CTA observations
        if self._method != 'NODES':
            if n_unbinned == 0 and n_binned != 0 and n_onoff == 0:
                self._binned_mode = True
                self._method      = 'SLICE'
            elif n_unbinned == 0 and n_binned == 0 and n_onoff != 0:
                self._onoff_mode = True
                self._method      = 'SLICE'
            elif n_unbinned == 0 and n_binned != 0 and n_onoff != 0:
                msg = 'Mix of binned and On/Off CTA observations found in ' \
                      'observation container. csscript does not support this ' \
                      'mix.'
                raise RuntimeError(msg)
            elif n_unbinned != 0 and (n_binned != 0 or n_onoff != 0):
                msg = 'Mix of unbinned and binned or On/Off CTA observations ' \
                      'found in observation container. csscript does not support ' \
                      'this mix.'
                raise RuntimeError(msg)
            elif n_unbinned != 0:
                self._method = 'SLICE'

        # Set ebounds
        self._set_ebounds()

        # Query other parameeters
        self['edisp'].boolean()
        self['calc_ulim'].boolean()
        self['calc_ts'].boolean()
        self['fix_bkg'].boolean()
        self['fix_srcs'].boolean()

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Write spectrum method header and parameters
        self._log_header1(gammalib.TERSE, 'Spectrum method')
        self._log_value(gammalib.TERSE, 'Unbinned CTA observations', n_unbinned)
        self._log_value(gammalib.TERSE, 'Binned CTA observations', n_binned)
        self._log_value(gammalib.TERSE, 'On/off CTA observations', n_onoff)
        self._log_value(gammalib.TERSE, 'Other observations', n_other)

        # If there are a mix of CTA and non-CTA observations and the method
        # is 'SLICE' then log a warning that non-CTA observations will be
        # ignored
        if n_cta > 0 and n_other > 0 and self._method == 'SLICE':
            self._log_string(gammalib.TERSE, '\nWARNING: Only CTA observation '
                             'can be handled with the "SLICE" method, all '
                             'non-CTA observation will be ignored.')

        # If there are only non-CTA observations and the method is 'SLICE'
        # then stop now
        elif n_other > 0:
            if self._method == 'SLICE':
                msg = 'Selected "SLICE" method but none of the observations ' \
                      'is a CTA observation. Please select "AUTO" or "NODES" ' \
                      'if no CTA observation is provided.'
                raise RuntimeError(msg)
            else:
                self._method = 'NODES'

        # Log selected spectrum method
        self._log_value(gammalib.TERSE, 'Selected spectrum method', self._method)

        # Return
        return

    def _set_ebounds(self):
        """
        Set energy boundaries
        """
        # If we are in binned or On/Off mode then align the energy boundaries
        # with the cube of RMF spectrum
        if self._binned_mode or self._onoff_mode:

            # Initialise energy boundaries for spectrum
            self._ebounds = gammalib.GEbounds()

            # Create energy boundaries according to user parameters
            ebounds = self._create_ebounds()

            # Extract binned energy boundaries from first observation in
            # container. This assumes that all observations have the same
            # energy binning. But even if this is not the case, the script
            # should work reasonably well since for each observation the
            # relevant energy bins will be selected.
            if self._binned_mode:
                cube_ebounds = self.obs()[0].events().ebounds()
            else:
                cube_ebounds = self.obs()[0].rmf().emeasured()

            # Loop over user energy boundaries and collect all energy bins
            # that overlap
            for i in range(ebounds.size()):

                # Extract minimum and maximum energy of user energy bin,
                # including some rounding tolerance
                emin = ebounds.emin(i).TeV() * 0.999 # Rounding tolerance
                emax = ebounds.emax(i).TeV() * 1.001 # Rounding tolerance

                # Set number of overlapping energy bins
                nbins = 0

                # Search first cube bin that is comprised within user energy
                # bin
                emin_value = -1.0
                for k in range(cube_ebounds.size()):
                    if cube_ebounds.emin(k).TeV() >= emin and \
                       cube_ebounds.emax(k).TeV() <= emax:
                        emin_value = cube_ebounds.emin(k).TeV()
                        break
                if emin_value < 0.0:
                    continue

                # Search last cube bin that is comprised within user energy bin
                for k in range(cube_ebounds.size()):
                    if cube_ebounds.emin(k).TeV() >= emin and \
                       cube_ebounds.emax(k).TeV() <= emax:
                        emax_value = cube_ebounds.emax(k).TeV()
                        nbins += 1

                # Append energy bin if there are overlapping bins in the
                # counts cube
                if nbins > 0:
                    self._ebounds.append(gammalib.GEnergy(emin_value, 'TeV'),
                                         gammalib.GEnergy(emax_value, 'TeV'))

            # Raise an exception if there are no overlapping layers
            if (len(self._ebounds) == 0):
                msg = 'Energy range ['+str(cube_ebounds.emin())+ \
                      ', '+str(cube_ebounds.emax())+'] of counts '+ \
                      'cube does not overlap with specified energy '+ \
                      'range ['+ \
                      str(ebounds.emin())+', '+str(ebounds.emax())+'].'+ \
                      ' Specify overlapping energy range.'
                raise RuntimeError(msg)

        # Unbinned mode
        else:
            self._ebounds = self._create_ebounds()

        # Return
        return

    def _log_spectral_binning(self):
        """
        Log spectral binning
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Spectral binning')

        # Log counts cube energy range for binned mode
        if self._binned_mode:
            cube_ebounds = self.obs()[0].events().ebounds()
            value = '%s - %s' % (str(cube_ebounds.emin()),
                                 str(cube_ebounds.emax()))
            self._log_value(gammalib.TERSE, 'Counts cube energy range', value)

        # Log RMF energy range for On/Off mode
        elif self._onoff_mode:
            etrue = self.obs()[0].rmf().etrue()
            ereco = self.obs()[0].rmf().emeasured()
            vtrue = '%s - %s' % (str(etrue.emin()), str(etrue.emax()))
            vreco = '%s - %s' % (str(ereco.emin()), str(ereco.emax()))
            self._log_value(gammalib.TERSE, 'RMF true energy range', vtrue)
            self._log_value(gammalib.TERSE, 'RMF measured energy range', vreco)

        # Log energy bins
        for i in range(self._ebounds.size()):
            value = '%s - %s' % (str(self._ebounds.emin(i)),
                                 str(self._ebounds.emax(i)))
            self._log_value(gammalib.TERSE, 'Bin %d' % (i+1), value)

        # Return
        return

    def _adjust_model_pars(self):
        """
        Adjust model parameters
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Adjust model parameters')

        # Adjust model parameters dependent on input user parameters
        for model in self.obs().models():

            # Initialise TS flag for all models to false
            model.tscalc(False)

            # Log model name
            self._log_header3(gammalib.EXPLICIT, model.name())

            # Deal with the source of interest    
            if model.name() == self['srcname'].string():

                # Fix all model parameters
                for par in model:
                    if par.is_free():
                        self._log_string(gammalib.EXPLICIT,
                                         ' Fixing "'+par.name()+'"')
                    par.fix()

                # Free the normalisation parameter which is assumed to be
                # the first spectral parameter
                normpar = model.spectral()[0]
                if normpar.is_fixed():
                    self._log_string(gammalib.EXPLICIT,
                                     ' Freeing "'+normpar.name()+'"')
                normpar.free()

                # Optionally compute Test Statistic value
                if self['calc_ts'].boolean():
                    model.tscalc(True)

            # Deal with background models
            elif self['fix_bkg'].boolean() and \
                 not model.classname() == 'GModelSky':
                for par in model:
                    if par.is_free():
                        self._log_string(gammalib.EXPLICIT,
                                         ' Fixing "'+par.name()+'"')
                    par.fix()

            # Deal with source models
            elif self['fix_srcs'].boolean() and \
                 model.classname() == 'GModelSky':
                for par in model:
                    if par.is_free():
                        self._log_string(gammalib.EXPLICIT,
                                         ' Fixing "'+par.name()+'"')
                    par.fix()

        # Return
        return

    def _set_replace_src_spectrum_by_nodes(self):
        """
        Replace source spectrum by node function
        """
        # Initialise model container
        models = gammalib.GModels()

        # Loop over model containers
        for model in self.obs().models():

            # If we deal with source model then replace the spectral model
            # by a node function
            if model.name() == self['srcname'].string():

                # Setup energies at log mean energies of bins
                energies = gammalib.GEnergies()
                for i in range(self._ebounds.size()):
                    energies.append(self._ebounds.elogmean(i))

                # Setup spectral node function
                spectrum = gammalib.GModelSpectralNodes(model.spectral(), energies)
                spectrum.autoscale()

                # Make sure that all nodes are positive
                for i in range(spectrum.nodes()):
                    par     = spectrum[i*2+1]
                    value   = par.value()
                    minimum = 1.0e-30 * value
                    if minimum <= 0.0:
                        minimum = 1.0e-50
                        if minimum < value:
                            value = minimum
                    par.value(value)
                    par.min(minimum)

                # Set spectral component of source model
                model.spectral(spectrum)

                # Append model
                models.append(model)

            # ... otherwise just append model
            else:
                models.append(model)

        # Put new model in observation containers
        self.obs().models(models)

        # Return
        return

    def _select_onoff_obs(self, obs, emin, emax):
        """
        Select an energy interval from one CTA On/Off observation

        Parameters
        ----------
        obs : `~gammalib.GCTAOnOffObservation`
            Minimum energy
        emin : `~gammalib.GEnergy()`
            Minimum energy
        emax : `~gammalib.GEnergy()`
            Maximum energy

        Returns
        -------
        obs : `~gammalib.GCTAOnOffObservation`
            CTA On/Off observation
        """
        # Select energy bins in etrue and ereco. 20% margin is added for
        # true energies to take into account the energy dispersion. 0.1%
        # margin is added for reconstructed energies to accomodate for
        # rounding errors.
        etrue = gammalib.GEbounds()
        ereco = gammalib.GEbounds()
        itrue = []
        ireco = []
        for i in range(obs.rmf().etrue().size()):
            etrue_min = obs.rmf().etrue().emin(i)
            etrue_max = obs.rmf().etrue().emax(i)
            if etrue_min * 1.2 >= emin and etrue_max * 0.8 <= emax:
                etrue.append(etrue_min, etrue_max)
                itrue.append(i)
        for i in range(obs.rmf().emeasured().size()):
            ereco_min = obs.rmf().emeasured().emin(i)
            ereco_max = obs.rmf().emeasured().emax(i)
            if ereco_min * 1.001 >= emin and ereco_max * 0.999 <= emax:
                ereco.append(ereco_min, ereco_max)
                ireco.append(i)

        # Extract PHA
        pha_on  = gammalib.GPha(ereco)
        pha_off = gammalib.GPha(ereco)
        pha_on.exposure(obs.on_spec().exposure())
        pha_off.exposure(obs.on_spec().exposure())
        for idst, isrc in enumerate(ireco):
            # On
            pha_on[idst] = obs.on_spec()[isrc]
            pha_on.areascal(idst, obs.on_spec().areascal(isrc))
            pha_on.backscal(idst, obs.on_spec().backscal(isrc))
            # Off
            pha_off[idst] = obs.off_spec()[isrc]
            pha_off.areascal(idst, obs.off_spec().areascal(isrc))
            pha_off.backscal(idst, obs.off_spec().backscal(isrc))

        # Extract ARF
        arf          = gammalib.GArf(etrue)
        arf_backresp = obs.arf()['BACKRESP']
        backresp     = []
        for idst, isrc in enumerate(itrue):
            arf[idst] = obs.arf()[isrc]
            backresp.append(arf_backresp[isrc])
        arf.append('BACKRESP', backresp)

        # Extract RMF
        rmf = gammalib.GRmf(etrue, ereco)
        for idst_true, isrc_true in enumerate(itrue):
            for idst_reco, isrc_reco in enumerate(ireco):
                rmf[idst_true, idst_reco] = obs.rmf()[isrc_true, isrc_reco]

        # Set On/Off observations
        obsid     = obs.id()
        statistic = obs.statistic()
        obs = gammalib.GCTAOnOffObservation(pha_on, pha_off, arf, rmf)
        obs.id(obsid)
        obs.statistic(statistic)

        # Return observation
        return obs

    def _select_obs(self, emin, emax):
        """
        Select observations for energy interval

        Parameters
        ----------
        emin : `~gammalib.GEnergy()`
            Minimum energy
        emax : `~gammalib.GEnergy()`
            Maximum energy

        Returns
        -------
        obs : `~gammalib.GObservations`
            Observation container
        """
        # Use ctcubemask for binned analysis
        if self._binned_mode:

            # Write header
            self._log_header3(gammalib.EXPLICIT, 'Filtering cube')

            # Select layers
            cubemask            = ctools.ctcubemask(self.obs())
            cubemask['regfile'] = 'NONE'
            cubemask['ra']      = 'UNDEFINED'
            cubemask['dec']     = 'UNDEFINED'
            cubemask['rad']     = 'UNDEFINED'
            cubemask['emin']    = emin.TeV()
            cubemask['emax']    = emax.TeV()
            cubemask.run() 

            # Set new binned observation
            obs = cubemask.obs().copy()

        # Use ...
        elif self._onoff_mode:

            # Write header
            self._log_header3(gammalib.EXPLICIT, 'Filtering PHA, ARF and RMF')

            # Initialise observation container
            obs = gammalib.GObservations()

            # Loop over all input observations and select energy bins for
            # all On/Off observations
            for run in self.obs():
                if run.classname() == 'GCTAOnOffObservation':
                    obs.append(self._select_onoff_obs(run, emin, emax))

            # Append models
            obs.models(self.obs().models())

        # Use ctselect for unbinned analysis
        else:

            # Write header
            self._log_header3(gammalib.EXPLICIT, 'Selecting events')

            # Select events
            select = ctools.ctselect(self.obs())
            select['ra']   = 'UNDEFINED'
            select['dec']  = 'UNDEFINED'
            select['rad']  = 'UNDEFINED'
            select['emin'] = emin.TeV()
            select['emax'] = emax.TeV()
            select['tmin'] = 'UNDEFINED'
            select['tmax'] = 'UNDEFINED'
            select.run()

            # Retrieve observation
            obs = select.obs().copy()

        # Return observation container
        return obs

    def _fit_energy_bin(self, i):
        """
        Fit data for one energy bin

        Parameters
        ----------
        i : int
            Energy bin index

        Returns
        -------
        result : dict
            Dictionary with fit results
        """
        # Get energy boundaries
        emin      = self._ebounds.emin(i)
        emax      = self._ebounds.emax(i)
        elogmean  = self._ebounds.elogmean(i)

        # Select observations for energy bin
        obs = self._select_obs(emin, emax)

        # Initialise dictionary
        result = {'energy':      elogmean.TeV(),
                  'energy_low':  (elogmean - emin).TeV(),
                  'energy_high': (emax - elogmean).TeV(),
                  'flux':        0.0,
                  'flux_err':    0.0,
                  'TS':          0.0,
                  'ulimit':      0.0,
                  'Npred':       0.0}

        # Write header for fitting
        self._log_header3(gammalib.EXPLICIT, 'Performing fit')

        # Perform maximum likelihood fit
        like          = ctools.ctlike(obs)
        like['edisp'] = self['edisp'].boolean()
        like.run()

        # Continue only if log-likelihood is non-zero
        if like.obs().logL() != 0.0:

            # Get results
            fitted_models = like.obs().models()
            source        = fitted_models[self['srcname'].string()]

            # Extract Test Statistic value
            if self['calc_ts'].boolean():
                result['TS'] = source.ts()

            # Compute Npred value (only works for unbinned analysis)
            if not self._binned_mode and not self._onoff_mode:
                for observation in like.obs():
                    result['Npred'] += observation.npred(source)

            # Compute upper flux limit
            ulimit_value = -1.0
            if self['calc_ulim'].boolean():

                # Logging information
                self._log_header3(gammalib.EXPLICIT, 'Computing upper limit')

                # Create upper limit object  
                ulimit = ctools.ctulimit(like.obs())
                ulimit['srcname'] = self['srcname'].string()
                ulimit['eref']    = elogmean.TeV()

                # Try to run upper limit and catch exceptions
                try:
                    ulimit.run()
                    ulimit_value = ulimit.diff_ulimit()
                except:
                    self._log_string(gammalib.EXPLICIT, 'Upper limit '
                                     'calculation failed.')
                    ulimit_value = -1.0

                # Compute upper limit
                if ulimit_value > 0.0:
                    result['ulimit'] = ulimit_value * elogmean.MeV() * \
                                       elogmean.MeV() * gammalib.MeV2erg

            # Compute differential flux and flux error
            fitted_flux = source.spectral().eval(elogmean)
            parvalue    = source.spectral()[0].value()
            if parvalue != 0.0:
                rel_error = source.spectral()[0].error() / parvalue
                e_flux    = fitted_flux * rel_error
            else:
                e_flux = 0.0

            # Convert differential flux and flux error to nuFnu
            elogmean2          = elogmean.MeV() * elogmean.MeV()
            result['flux']     = fitted_flux * elogmean2 * gammalib.MeV2erg
            result['flux_err'] = e_flux      * elogmean2 * gammalib.MeV2erg

            # Log information
            value = '%e +/- %e' % (result['flux'], result['flux_err'])
            if self['calc_ulim'].boolean() and result['ulimit'] > 0.0:
                value += ' [< %e]' % (result['ulimit'])
            value += ' erg/cm2/s'
            if self['calc_ts'].boolean() and result['TS'] > 0.0:
                value += ' (TS = %.3f)' % (result['TS'])
            self._log_value(gammalib.TERSE, 'Bin '+str(i+1), value)

        # ... otherwise if logL is zero then signal that bin is
        # skipped
        else:
            value = 'No event in this bin. Likelihood is zero. Bin is skipped.'
            self._log_value(gammalib.TERSE, 'Bin '+str(i+1), value)

        # Return result
        return result

    def _fit_energy_bins(self):
        """
        Fit model to energy bins

        Returns
        -------
        results : list of dict
            List of dictionaries with fit results
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Generate spectrum')
        self._log_string(gammalib.TERSE, str(self._ebounds))

        # Initialise results
        results = []

        # Loop over energy bins
        for i in range(self._ebounds.size()):

            # Write header for energy bin
            self._log_header2(gammalib.EXPLICIT, 'Energy bin '+str(i+1))

            # Fit energy bin
            result = self._fit_energy_bin(i)

            # Append results
            results.append(result)

        # Return results
        return results

    def _fit_model(self):
        """
        Fit model to observations

        Returns
        -------
        results : list of dict
            List of dictionaries with fit results
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Generate spectrum')

        # Write header for fitting
        self._log_header3(gammalib.EXPLICIT, 'Performing fit')

        # Perform maximum likelihood fit
        like          = ctools.ctlike(self.obs())
        like['edisp'] = self['edisp'].boolean()
        like.run()

        # Initialise fit results
        results = []

        # Extract fit results
        model    = like.obs().models()[self['srcname'].string()]
        spectrum = model.spectral()
        ts       = model.ts()

        # Loop over all nodes
        for i in range(spectrum.nodes()):

            # Get energy boundaries
            emin      = self._ebounds.emin(i)
            emax      = self._ebounds.emax(i)
            elogmean  = self._ebounds.elogmean(i)

            # Initialise dictionary
            result = {'energy':      elogmean.TeV(),
                      'energy_low':  (elogmean - emin).TeV(),
                      'energy_high': (emax - elogmean).TeV(),
                      'flux':        0.0,
                      'flux_err':    0.0,
                      'TS':          0.0,
                      'ulimit':      0.0,
                      'Npred':       0.0}

            # Convert differential flux and flux error to nuFnu
            norm               = elogmean.MeV() * elogmean.MeV()  * gammalib.MeV2erg
            result['flux']     = spectrum[i*2+1].value() * norm
            result['flux_err'] = spectrum[i*2+1].error() * norm

            # Compute TS
            if self['calc_ts'].boolean():

                # Copy observation container
                obs = like.obs().copy()

                # Set intensity of node to tiny value
                par = obs.models()[self['srcname'].string()].spectral()[i*2+1]
                value = 1.0e-30 * par.value()
                if par.min() > value:
                    par.min(value)
                par.value(value)
                par.fix()

                # Perform maximum likelihood fit
                tslike          = ctools.ctlike(obs)
                tslike['edisp'] = self['edisp'].boolean()
                tslike.run()

                # Store Test Statistic
                model        = tslike.obs().models()[self['srcname'].string()]
                result['TS'] = ts - model.ts()

            # Log information
            value = '%e +/- %e' % (result['flux'], result['flux_err'])
            if self['calc_ulim'].boolean() and result['ulimit'] > 0.0:
                value += ' [< %e]' % (result['ulimit'])
            value += ' erg/cm2/s'
            if self['calc_ts'].boolean() and result['TS'] > 0.0:
                value += ' (TS = %.3f)' % (result['TS'])
            self._log_value(gammalib.TERSE, 'Bin '+str(i+1), value)

            # Append results
            results.append(result)

        # Return results
        return results

    def _create_fits(self, results):
        """
        Create FITS file

        Parameters
        ----------
        results : list of dict
            List of result dictionaries
        """
        # Create FITS table columns
        nrows        = self._ebounds.size()
        energy       = gammalib.GFitsTableDoubleCol('Energy', nrows)
        energy_low   = gammalib.GFitsTableDoubleCol('ed_Energy', nrows)
        energy_high  = gammalib.GFitsTableDoubleCol('eu_Energy', nrows)
        flux         = gammalib.GFitsTableDoubleCol('Flux', nrows)
        flux_err     = gammalib.GFitsTableDoubleCol('e_Flux', nrows)
        TSvalues     = gammalib.GFitsTableDoubleCol('TS', nrows)
        ulim_values  = gammalib.GFitsTableDoubleCol('UpperLimit', nrows)
        Npred_values = gammalib.GFitsTableDoubleCol('Npred', nrows)
        energy.unit('TeV')
        energy_low.unit('TeV')
        energy_high.unit('TeV')
        flux.unit('erg/cm2/s')
        flux_err.unit('erg/cm2/s')
        ulim_values.unit('erg/cm2/s')

        # File FITS table columns
        for i, result in enumerate(results):
            energy[i]       = result['energy']
            energy_low[i]   = result['energy_low']
            energy_high[i]  = result['energy_high']
            flux[i]         = result['flux']
            flux_err[i]     = result['flux_err']
            TSvalues[i]     = result['TS']
            ulim_values[i]  = result['ulimit']
            Npred_values[i] = result['Npred']

        # Initialise FITS Table with extension "SPECTRUM"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('SPECTRUM')

        # Add Header for compatibility with gammalib.GMWLSpectrum
        table.card('INSTRUME', 'CTA', 'Name of Instrument')
        table.card('TELESCOP', 'CTA', 'Name of Telescope')

        # Append filled columns to fits table    
        table.append(energy)
        table.append(energy_low)
        table.append(energy_high)
        table.append(flux)
        table.append(flux_err)
        table.append(TSvalues)
        table.append(ulim_values)
        table.append(Npred_values)

        # Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append(table)

        # Return
        return


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

        # Write input observation container into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Input observation')

        # Write spectral binning into logger
        self._log_spectral_binning()

        # Adjust model parameters dependent on input user parameters
        self._adjust_model_pars()

        # Case A: "SLICE" method
        if self._method == 'SLICE':

            # Fit energy bins
            results = self._fit_energy_bins()

        # Case B: "NODES" method
        elif self._method == 'NODES':

            # Replace source spectrum by nodes function
            self._set_replace_src_spectrum_by_nodes()

            # Fit model
            results = self._fit_model()

        # Create FITS file
        self._create_fits(results)

        # Optionally publish spectrum
        if self['publish'].boolean():
            self.publish()

        # Return
        return

    def save(self):
        """
        Save spectrum
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save spectrum')

        # Continue only if FITS file is valid
        if self._fits != None:

            # Get outmap parameter
            outfile = self['outfile'].filename()

            # Log file name
            self._log_value(gammalib.NORMAL, 'Spectrum file', outfile.url())

            # Save spectrum
            self._fits.saveto(outfile, self['clobber'].boolean())

        # Return
        return

    def publish(self, name=''):
        """
        Publish spectrum

        Parameters
        ----------
        name : str, optional
            Name of spectrum
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Publish spectrum')

        # Continue only if FITS file is valid
        if self._fits != None:

            # Continue only if spectrum is valid
            if self._fits.contains('SPECTRUM'):

                # Set default name is user name is empty
                if not name:
                    user_name = self._name()
                else:
                    user_name = name

                # Log file name
                self._log_value(gammalib.NORMAL, 'Spectrum name', user_name)

                # Publish spectrum
                self._fits.publish('SPECTRUM', user_name)

        # Return
        return

    def spectrum(self):
        """
        Return spectrum FITS file

        Returns:
            FITS file containing spectrum
        """
        # Return
        return self._fits

    def models(self, models):
        """
        Set model
        """
        # Copy models
        self.obs().models(models.clone())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csspec(sys.argv)

    # Execute application
    app.execute()
