#! /usr/bin/env python
# ==========================================================================
# Computes the PHA spectra for source/background and ARF/RMF files using the
# reflected region method
#
# Copyright (C) 2017-2018 Luigi Tibaldo
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
import math
import sys
from cscripts import mputils


# =============== #
# csfindobs class #
# =============== #
class csphagen(ctools.csobservation):
    """
    Generate PHA, ARF and RMF files for classical IACT spectral analysis
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise other variables
        self._ebounds       = gammalib.GEbounds()
        self._etruebounds   = gammalib.GEbounds()
        self._src_dir       = gammalib.GSkyDir()
        self._src_reg       = gammalib.GSkyRegions()
        self._models        = gammalib.GModels()
        self._srcname       = ''
        self._use_model_bkg = True
        self._bkg_regs      = []
        self._excl_reg      = None
        self._has_exclusion = False
        self._srcshape      = ''
        self._rad           = 0.0
        self._nthreads      = 0

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
        state = {'base'          : ctools.csobservation.__getstate__(self),
                 'ebounds'       : self._ebounds,
                 'etruebounds'   : self._etruebounds,
                 'src_dir'       : self._src_dir,
                 'src_reg'       : self._src_reg,
                 'models'        : self._models,
                 'srcname'       : self._srcname,
                 'use_model_bkg' : self._use_model_bkg,
                 'bkg_regs'      : self._bkg_regs,
                 'excl_reg'      : self._excl_reg,
                 'has_exclusion' : self._has_exclusion,
                 'srcshape'      : self._srcshape,
                 'rad'           : self._rad,
                 'nthreads'      : self._nthreads}

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
        self._ebounds       = state['ebounds']
        self._etruebounds   = state['etruebounds']
        self._src_dir       = state['src_dir']
        self._src_reg       = state['src_reg']
        self._models        = state['models']
        self._srcname       = state['srcname']
        self._use_model_bkg = state['use_model_bkg']
        self._bkg_regs      = state['bkg_regs']
        self._excl_reg      = state['excl_reg']
        self._has_exclusion = state['has_exclusion']
        self._srcshape      = state['srcshape']
        self._rad           = state['rad']
        self._nthreads      = state['nthreads']

        # Return
        return

    # Private methods
    def _query_src_direction(self):
        """
        Set up the source direction parameter
        """
        # Initialise source direction
        self._src_dir = gammalib.GSkyDir()

        # Get coordinate systel
        coordsys = self['coordsys'].string()

        # If coordinate system is celestial then query "ra" and "dec"
        if coordsys == 'CEL':
            ra  = self['ra'].real()
            dec = self['dec'].real()
            self._src_dir.radec_deg(ra, dec)

        # ... otherwise, if coordinate system is galactic then query "glon"
        # and "glat"
        elif coordsys == 'GAL':
            glon = self['glon'].real()
            glat = self['glat'].real()
            self._src_dir.lb_deg(glon, glat)

        # Return
        return

    def _get_parameters_bkgmethod_reflected(self):
        """
        Get parameters for REFLECTED background method
        """
        # Get source shape
        self._srcshape = self['srcshape'].string()

        # Set source region (so far only CIRCLE is supported)
        if self._srcshape == 'CIRCLE':

            # Query source direction
            self._query_src_direction()

            # Query minimum number of background regions
            self['bkgregmin'].integer()

            # Set circular source region
            self._rad = self['rad'].real()
            self._src_reg.append(gammalib.GSkyRegionCircle(self._src_dir, self._rad))

        # Query usage of background model
        self['use_model_bkg'].boolean()

        # Return
        return

    def _get_parameters_bkgmethod_custom(self):
        """
        Get parameters for CUSTOM background method

        Raises
        ------
        RuntimeError
            Only one On region is allowed
        """
        # Set up source region
        filename      = self['srcregfile'].filename()
        self._src_reg = gammalib.GSkyRegions(filename)

        # Raise an exception if there is more than one source region
        if len(self._src_reg) != 1:
            raise RuntimeError('Only one On region is allowed')

        # Set up source direction. Query parameters if neccessary.
        if isinstance(self._src_reg[0], gammalib.GSkyRegionCircle):
            self._src_dir = self._src_reg[0].centre()
            self._rad     = self._src_reg[0].radius()
        else:
            self._query_src_direction()

        # Make sure that all CTA observations have an Off region by loading the
        # Off region region the parameter 'bkgregfile' for all CTA observations
        # without Off region
        for obs in self.obs():
            if obs.classname() == 'GCTAObservation':
                if obs.off_regions().is_empty():
                    filename = self['bkgregfile'].filename()
                    regions  = gammalib.GSkyRegions(filename)
                    obs.off_regions(regions)

        # Return
        return

    def _get_parameters_bkgmethod(self):
        """
        Get background method parameters
        """
        # Get background method
        bkgmethod = self['bkgmethod'].string()

        # Get background method dependent parameters
        if bkgmethod == 'REFLECTED':
            self._get_parameters_bkgmethod_reflected()
        elif bkgmethod == 'CUSTOM':
            self._get_parameters_bkgmethod_custom()

        # Query parameters that are needed for all background methods
        self['maxoffset'].real()

        # Return
        return

    def _get_parameters(self):
        """
        Get parameters from parfile and setup observations
        """
        # Clear source models
        self._models.clear()

        # Setup observations (require response and allow event list, don't
        # allow counts cube)
        self._setup_observations(self.obs(), True, True, False)

        # Get source model and source name. First try to extract models from
        # observation container. If this does not work then try creating
        # model from the inmodel parameter
        if self.obs().models().size() > 0:
            self._models  = self.obs().models().clone()
            self._srcname = self['srcname'].string()
        elif self['inmodel'].is_valid():
            inmodel         = self['inmodel'].filename()
            self._models    = gammalib.GModels(inmodel)
            self._srcname   = self['srcname'].string()

        # Set energy bounds
        self._ebounds = self._create_ebounds()

        # Initialize empty src regions container
        self._src_reg = gammalib.GSkyRegions()

        # Exclusion map
        if self['inexclusion'].is_valid():
            inexclusion         = self['inexclusion'].filename()
            self._excl_reg      = gammalib.GSkyRegionMap(inexclusion)
            self._has_exclusion = True
        else:
            self._has_exclusion = False

        # Query remaining parametersStacking
        self['stack'].boolean()

        # Query ahead output parameters
        if (self._read_ahead()):
            self['outobs'].string()
            self['prefix'].string()

        # If there are no observations in container then get them from the
        # parameter file
        if self.obs().is_empty():
            self.obs(self._get_observations(False))

        # Get background method parameters (have to come after setting up of
        # observations)
        self._get_parameters_bkgmethod()

        # Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Set number of processes for multiprocessing
        self._nthreads = mputils.nthreads(self)

        # If we have no model then create now a dummy model
        if self._models.is_empty():
            spatial  = gammalib.GModelSpatialPointSource(self._src_dir)
            spectral = gammalib.GModelSpectralPlaw()
            model    = gammalib.GModelSky(spatial, spectral)
            model.name('Dummy')
            self._models.append(model)
            self._srcname       = 'Dummy'
            self._use_model_bkg = False

        # Return
        return

    def _reflected_regions(self, obs):
        """
        Calculate list of reflected regions for a single observation (pointing)

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation()`
            CTA observation

        Returns
        -------
        regions : `~gammalib.GSkyRegions`
            List of reflected regions
        """
        # Initialise list of reflected regions
        regions = gammalib.GSkyRegions()

        # Get offset angle of source
        pnt_dir = obs.pointing().dir()
        offset  = pnt_dir.dist_deg(self._src_dir)

        # If ...
        if offset <= self._rad or offset >= self['maxoffset'].real():
            self._log_string(gammalib.EXPLICIT, 'Observation %s pointed at %.3f '
                             'deg from source' % ((obs.id(), offset)))

        # ... otherwise
        else:
            posang = pnt_dir.posang_deg(self._src_dir)
            if self._srcshape == 'CIRCLE':
                # angular separation of reflected regions wrt camera center
                # and number
                alpha = 1.05 * 2.0 * self._rad / offset
                # 1.05 ensures background regions do not overlap due to
                # numerical precision issues
                N = int(2.0 * math.pi / alpha)
                if N < self['bkgregmin'].integer() + 3:
                    self._log_string(gammalib.EXPLICIT, 'Observation %s: '
                                     'insufficient regions for background '
                                     'estimation' % (obs.id()))

                # Otherwise loop over position angle to create reflected
                # regions
                else:
                    alpha = 360.0 / N
                    for s in range(2, N - 1):
                        dphi    = s * alpha
                        ctr_dir = pnt_dir.clone()
                        ctr_dir.rotate_deg(posang + dphi, offset)
                        region = gammalib.GSkyRegionCircle(ctr_dir, self._rad)
                        if self._has_exclusion:
                            if self._excl_reg.overlaps(region):
                                self._log_string(gammalib.VERBOSE, 'Observation '
                                                 '%s: reflected region overlaps '
                                                 'with exclusion region' %
                                                 (obs.id()))
                            else:
                                regions.append(region)
                        else:
                            regions.append(region)

        # Return reflected regions
        return regions

    def _set_models(self, obs):
        """
        Set models in observation container

        The method replaces all "CTA", "HESS", "VERITAS", and "MAGIC"
        background models by "CTAOnOff" background models.

        Parameters
        ----------
        obs : `~gammalib.GObservations()`
            Observation container

        Returns
        -------
        obs : `~gammalib.GObservations()`
            Observation container
        """
        # Initialise model container
        models = gammalib.GModels()

        # Loop over all models and replace CTA/HESS/VERITAS/MAGIC background
        # models by CTAOnOff background model
        for model in self.obs().models():
            if 'GCTA' in model.classname():
                if 'CTA'     in model.instruments() or \
                   'HESS'    in model.instruments() or \
                   'VERITAS' in model.instruments() or \
                   'MAGIC'   in model.instruments():
                    model.instruments('CTAOnOff')
            models.append(model)

        # Append model to observation container
        obs.models(models)

        # Return observation container
        return obs

    def _etrue_ebounds(self):
        """
        Set true energy boundaries

        Returns
        -------
        ebounds : `~gammalib.GEbounds()`
            True energy boundaries
        """
        # Determine minimum and maximum energies
        emin = self._ebounds.emin() * 0.5
        emax = self._ebounds.emax() * 1.2
        if emin.TeV() < self['etruemin'].real():
            emin = gammalib.GEnergy(self['etruemin'].real(), 'TeV')
        if emax.TeV() > self['etruemax'].real():
            emax = gammalib.GEnergy(self['etruemax'].real(), 'TeV')

        # Determine number of energy bins
        n_decades = (emax.log10TeV() - emin.log10TeV())
        n_bins    = int(n_decades * float(self['etruebins'].integer()) + 0.5)
        if n_bins < 1:
            n_bins = 1

        # Set energy boundaries
        ebounds = gammalib.GEbounds(n_bins, emin, emax)

        # Write header
        self._log_header1(gammalib.TERSE, 'True energy binning')

        # Log true energy bins
        for i in range(ebounds.size()):
            value = '%s - %s' % (str(ebounds.emin(i)), str(ebounds.emax(i)))
            self._log_value(gammalib.TERSE, 'Bin %d' % (i+1), value)

        # Return energy boundaries
        return ebounds

    def _set_background_regions(self, obs):
        """
        Set background regions for an observation

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation()`
            CTA observation

        Returns
        -------
        regions : `~gammalib.GSkyRegions()`
            Background regions
        """
        # Initialise empty background regions for this observation
        bkg_reg = gammalib.GSkyRegions()

        # If reflected background is requested then created reflected
        # background regions
        if self['bkgmethod'].string() == 'REFLECTED':
            bkg_reg = self._reflected_regions(obs)

        # ... otherwise if custom background is requested then get the
        # background regions from the observation. We use a copy here since
        # otherwise the background regions go out of scope once the observations
        # are replaced by the On/Off observations.
        elif self['bkgmethod'].string() == 'CUSTOM':
            bkg_reg = obs.off_regions().copy()

        # Return background regions
        return bkg_reg

    def _process_observation(self,i):
        """
        Generate On/Off spectra for individual observation

        Parameters
        ----------
        i : int
            Observation number

        Returns
        -------
        result : dict
            On/Off spectra, background regions, observation id
        """
        # Retrieve observation from container
        onoff   = None
        bkg_reg = None
        obs     = self.obs()[i]

        # Skip non CTA observations
        if obs.classname() != 'GCTAObservation':
            self._log_string(gammalib.NORMAL, 'Skip %s observation "%s"' % \
                             (obs.instrument(), obs.id()))

        # Otherwise calculate On/Off spectra
        else:
            # Log current observation
            self._log_string(gammalib.NORMAL, 'Process observation "%s"' % \
                             (obs.id()))

            # Set background regions for this observation
            bkg_reg = self._set_background_regions(obs)

            # Get IRF background template usage flag
            use_model_bkg = self['use_model_bkg'].boolean()
            if not use_model_bkg:
                self._log_string(gammalib.NORMAL, 'Background model not used, '
                                 'assume constant backround rate.')

            # If there are background regions then create On/Off observation
            # and append it to the output container
            if bkg_reg.size() >= self['bkgregmin'].integer():
                onoff = gammalib.GCTAOnOffObservation(obs,
                                                      self._models,
                                                      self._srcname,
                                                      self._etruebounds,
                                                      self._ebounds,
                                                      self._src_reg,
                                                      bkg_reg,
                                                      use_model_bkg)
                onoff.id(obs.id())

            else:
                self._log_string(gammalib.NORMAL, 'Observation %s not included '
                                 'in spectra generation' % (obs.id()))

        # Construct dictionary with results
        result = {'onoff'   : onoff,
                  'bkg_reg' : bkg_reg,
                  'id'      : obs.id()}

        # Return results
        return result

    def _unpack_result(self, outobs, result):
        """
        Unpack result from calculation of On/Off regions

        Parameters
        ----------
        outobs : `~gammalib.GObservations`
            Observation container
        result : dict
            On/Off spectra, background regions, observation id

        Returns
        -------
        outobs : `~gammalib.GObservations`
            Observation container with result appended
        """
        # Continue only if result is valid
        if result['onoff'] != None:
        
            # If the results contain an On/Off observation
            if result['onoff'].classname() == 'GCTAOnOffObservation':

                # Append observation to observation container
                outobs.append(result['onoff'])

                # Append background regions
                self._bkg_regs.append({'regions': result['bkg_reg'],
                                      'id': result['id']})

        # Return observation container
        return outobs


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

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Observation')

        # Set true energy bins
        self._etruebounds = self._etrue_ebounds()

        # Write header
        self._log_header1(gammalib.TERSE, 'Spectral binning')

        # Log reconstructed energy bins
        for i in range(self._ebounds.size()):
            value = '%s - %s' % (str(self._ebounds.emin(i)),
                                 str(self._ebounds.emax(i)))
            self._log_value(gammalib.TERSE, 'Bin %d' % (i+1), value)

        # Write header
        self._log_header1(gammalib.NORMAL,
                          'Generation of source and background spectra')

        # Initialise run variables
        outobs         = gammalib.GObservations()
        self._bkg_regs = []

        # If there is more than one observation and we use multiprocessing
        if self._nthreads > 1 and self.obs().size() > 1:

            # Compute observations
            args        = [(self, '_process_observation', i)
                           for i in range(self.obs().size())]
            poolresults = mputils.process(self._nthreads, mputils.mpfunc, args)

            # Construct results
            for i in range(self.obs().size()):
                outobs = self._unpack_result(outobs,poolresults[i][0])
                self._log_string(gammalib.TERSE, poolresults[i][1]['log'], False)

        # Otherwise, loop through observations and generate pha, arf, rmf files
        else:
            for i in range(self.obs().size()):
                # Process individual observation
                result = self._process_observation(i)
                outobs = self._unpack_result(outobs, result)

        # Stack observations
        if outobs.size() > 1 and self['stack'].boolean() == True:

            # Write header
            self._log_header1(gammalib.NORMAL, 'Stacking %d observations' %
                              (outobs.size()))

            # Stack observations
            stacked_obs = gammalib.GCTAOnOffObservation(outobs)

            # Put stacked observations in output container
            outobs = gammalib.GObservations()
            outobs.append(stacked_obs)

        # Set models in output container
        outobs = self._set_models(outobs)

        # Set observation container
        self.obs(outobs)

        # Return
        return

    def save(self):
        """ 
        Save data
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save data')

        # Get XML output filename, prefix and clobber
        outobs  = self['outobs'].string()
        prefix  = self['prefix'].string()
        clobber = self['clobber'].boolean()

        # Loop over all observation in container
        for obs in self.obs():

            # Set filenames
            if self['stack'].boolean():
                onname  = prefix + '_stacked_pha_on.fits'
                offname = prefix + '_stacked_pha_off.fits'
                arfname = prefix + '_stacked_arf.fits'
                rmfname = prefix + '_stacked_rmf.fits'
            elif self.obs().size() > 1:
                onname  = prefix + '_%s_pha_on.fits' % (obs.id())
                offname = prefix + '_%s_pha_off.fits' % (obs.id())
                arfname = prefix + '_%s_arf.fits' % (obs.id())
                rmfname = prefix + '_%s_rmf.fits' % (obs.id())
            else:
                onname  = prefix + '_pha_on.fits'
                offname = prefix + '_pha_off.fits'
                arfname = prefix + '_arf.fits'
                rmfname = prefix + '_rmf.fits'

            # Set background and response file names in On spectrum
            obs.on_spec().backfile(offname)
            obs.on_spec().respfile(rmfname)
            obs.on_spec().ancrfile(arfname)

            # Save files
            obs.on_spec().save(onname, clobber)
            obs.off_spec().save(offname, clobber)
            obs.arf().save(arfname, clobber)
            obs.rmf().save(rmfname, clobber)

            # Log file names
            self._log_value(gammalib.NORMAL, 'PHA on file', onname)
            self._log_value(gammalib.NORMAL, 'PHA off file', offname)
            self._log_value(gammalib.NORMAL, 'ARF file', arfname)
            self._log_value(gammalib.NORMAL, 'RMF file', rmfname)

        # Save observation definition XML file
        self.obs().save(outobs)

        # Log file name
        self._log_value(gammalib.NORMAL, 'Obs. definition XML file', outobs)

        # Save ds9 On region file
        regname = prefix + '_on.reg'
        self._src_reg.save(regname)
        self._log_value(gammalib.NORMAL, 'On region file', regname)

        # Save ds9 Off region files
        for bkg_reg in self._bkg_regs:

            # Set filename
            if len(self._bkg_regs) > 1:
                regname = prefix + '_%s_off.reg' % (bkg_reg['id'])
            else:
                regname = prefix + '_off.reg'

            # Save ds9 region file
            bkg_reg['regions'].save(regname)

            # Log file name
            self._log_value(gammalib.NORMAL, 'Off region file', regname)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csphagen(sys.argv)

    # Execute application
    app.execute()
