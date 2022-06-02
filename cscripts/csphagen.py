#! /usr/bin/env python
# ==========================================================================
# Computes the PHA spectra for source/background and ARF/RMF files using the
# reflected region method
#
# Copyright (C) 2017-2022 Luigi Tibaldo
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
        self._obs_off       = gammalib.GObservations()
        self._ebounds       = gammalib.GEbounds()
        self._etruebounds   = gammalib.GEbounds()
        self._src_dir       = gammalib.GSkyDir()
        self._src_reg       = gammalib.GSkyRegions()
        self._models        = gammalib.GModels()
        self._srcname       = ''
        self._bkg_regs      = []
        self._excl_reg      = None
        self._has_exclusion = False
        self._srcshape      = ''
        self._rad           = 0.0
        self._reg_width     = 0.0
        self._reg_height    = 0.0
        self._reg_posang    = 0.0
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
                 'obs_off'       : self._obs_off,
                 'ebounds'       : self._ebounds,
                 'etruebounds'   : self._etruebounds,
                 'src_dir'       : self._src_dir,
                 'src_reg'       : self._src_reg,
                 'models'        : self._models,
                 'srcname'       : self._srcname,
                 'bkg_regs'      : self._bkg_regs,
                 'excl_reg'      : self._excl_reg,
                 'has_exclusion' : self._has_exclusion,
                 'srcshape'      : self._srcshape,
                 'rad'           : self._rad,
                 'reg_width'     : self._reg_width,
                 'reg_height'    : self._reg_height,
                 'reg_posang'    : self._reg_posang,
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
        self._obs_off       = state['obs_off']
        self._ebounds       = state['ebounds']
        self._etruebounds   = state['etruebounds']
        self._src_dir       = state['src_dir']
        self._src_reg       = state['src_reg']
        self._models        = state['models']
        self._srcname       = state['srcname']
        self._bkg_regs      = state['bkg_regs']
        self._excl_reg      = state['excl_reg']
        self._has_exclusion = state['has_exclusion']
        self._srcshape      = state['srcshape']
        self._rad           = state['rad']
        self._reg_width     = state['reg_width']
        self._reg_height    = state['reg_height']
        self._reg_posang    = state['reg_posang']
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

    def _compute_posang(self, pnt_dir, a, b):
        """
        Compute the difference in position angle wrt the pointing in degrees

        Parameters
        ----------
        pnt_dir : `~gammalib.GSkyDir`
            Pointing direction
        a : `~gammalib.GSkyDir`
            First sky direction
        a : `~gammalib.GSkyDir`
            Second sky direction

        Returns
        -------
        posang : float
            Position angle (degrees)
        """
        # Compute position angles
        posang_a = pnt_dir.posang_deg(a) % 360
        posang_b = pnt_dir.posang_deg(b) % 360

        # Compute difference
        posang = abs(posang_a - posang_b)

        # Return position angle
        return posang

    def _get_regions(self, filename):
        """
        Get regions from DS9 file or FITS file

        Parameters
        ----------
        filename : `~gammalib.GFilename`
            Filename

        Returns
        -------
        regs : `~gammalib.GSkyRegions`
            Region container
        """
        # If filename is a FITS file then load region map and append to
        # list of regions
        if filename.is_fits():
            map  = gammalib.GSkyRegionMap(filename)
            regs = gammalib.GSkyRegions()
            regs.append(map)

        # ... otherwise load DS9 file
        else:
            regs = gammalib.GSkyRegions(filename)

        # Return region container
        return regs

    def _get_source_parameters(self):
        """
        Get parameters to define source/On region
        """

        # Get source shape
        self._srcshape = self['srcshape'].string()

        # Query source direction
        self._query_src_direction()

        # If source shape is a circle the append GSkyRegionCircle
        if self._srcshape == 'CIRCLE':

            # Set circular source region
            self._rad = self['rad'].real()
            self._src_reg.append(gammalib.GSkyRegionCircle(self._src_dir, self._rad))

        # ... otherwise if source shape is a rectangle then append
        # GSkyRegionRectangle
        elif self._srcshape == 'RECT':

            # Set rectangular source region
            self._reg_width  = self['width'].real()
            self._reg_height = self['height'].real()
            self._reg_posang = self['posang'].real()
            self._src_reg.append(gammalib.GSkyRegionRectangle(self._src_dir,
                                                              self._reg_width,
                                                              self._reg_height,
                                                              self._reg_posang))

        # Return
        return

    def _get_parameters_bkgmethod_reflected(self):
        """
        Get parameters for REFLECTED background method
        """

        # Query parameters for source/On region definition
        self._get_source_parameters()

        # Query minimum number of background regions and
        # number of background regions to skip next to On region
        self['bkgregmin'].integer()
        self['bkgregskip'].integer()

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
        self._src_reg = self._get_regions(filename)

        # Raise an exception if there is more than one source region
        if len(self._src_reg) != 1:
            raise RuntimeError('Only one On region is allowed')

        # Set up source direction. Query parameters if neccessary.
        if self._models.is_empty():
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
                    regions  = self._get_regions(filename)
                    obs.off_regions(regions)

        # Return
        return

    def _get_parameters_bkgmethod_off(self):
        """
        Get parameters for OFF background method

        Raises
        ------
        RuntimeError
            On and Off observations must have same size
        RuntimeError
            Off observations must be event lists
        """

        # Set up Off observations. If there are no Off observations in the
        # container then load them via user parameters
        if self.obs_off().is_empty():

            # Get Off observation file name
            filename = self['inobsoff'].filename()

            # If Off observation is a FITS file then load observation and
            # append it to the Off observation container
            if gammalib.GFilename(filename).is_fits():
                self._obs_off.append(gammalib.GCTAObservation(filename))

            # ... otherwise load XML file into Off observation container
            else:
                self._obs_off.load(filename)

        # Check that size of On and Off observations are the same, otherwise
        # raise error
        if self.obs().size() != self.obs_off().size():
            raise RuntimeError('On and Off observations must have the same size')

        # Loop through observations
        for obs in self.obs_off():

            # Check that observation is event list, otherwise throw error
            if obs.eventtype() != "EventList":
                raise RuntimeError('Off observations must be event lists')

            # Check that they have response, otherwise assign based on user parameter
            if obs.has_response() == False:

                # Get database and IRF
                database = self["caldb"].string()
                irf      = self["irf"].string()

                # Create an XML element for response
                parameter = "parameter name=\"Calibration\"" +\
                               " database=\"" + database + "\"" +\
                               " response=\"" + irf + "\""
                xml = gammalib.GXmlElement()
                xml.append(parameter)

                # Create CTA response
                response = gammalib.GCTAResponseIrf(xml)

                # Attach response to observation
                obs.response(response)

        # Add models from Off observations to model container
        for model in self.obs_off().models():
            self._models.append(model)

        # Query parameters for source/On region definition
        self._get_source_parameters()

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
        elif bkgmethod == 'OFF':
            self._get_parameters_bkgmethod_off()

        # Query parameters that are needed for all background methods
        self['maxoffset'].real()
        self['use_model_bkg'].boolean()

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
            inmodel       = self['inmodel'].filename()
            self._models  = gammalib.GModels(inmodel)
            self._srcname = self['srcname'].string()

        # Set energy bounds
        self._ebounds = self._create_ebounds()

        # Initialize empty src regions container
        self._src_reg = gammalib.GSkyRegions()

        # Exclusion map
        if (self._excl_reg is not None) and (self._excl_reg.map().npix() > 0):
            # Exclusion map set and is not empty
            self._has_exclusion = True
        elif self['inexclusion'].is_valid():
            inexclusion         = self['inexclusion'].filename()
            # If the user has not specified the extension to use
            # and there is an extension called 'EXCLUSION' ...
            if not inexclusion.has_extname()\
                    and not inexclusion.has_extno()\
                    and gammalib.GFits(inexclusion).contains('EXCLUSION'):
                # ... choose it for the exclusion map
                extname = inexclusion.url() + '[EXCLUSION]'
                inexclusion = gammalib.GFilename(extname)
            # Otherwise will pick the default (primary) HDU
            self._excl_reg      = gammalib.GSkyRegionMap(inexclusion)
            self._has_exclusion = True
        else:
            self._has_exclusion = False

        # Get background method parameters (have to come after setting up of
        # observations and models)
        self._get_parameters_bkgmethod()

        # If there are multiple observations query whether to stack them
        if self.obs().size() > 1:
            self['stack'].boolean()

        # Query ahead output parameters
        if (self._read_ahead()):
            self['outobs'].filename()
            self['outmodel'].filename()
            self['prefix'].string()

        # Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Set number of processes for multiprocessing
        self._nthreads = mputils.nthreads(self)

        # If we have no model then create now a dummy model
        if self._models.is_empty():
            spatial  = gammalib.GModelSpatialPointSource(self._src_dir)
            spectral = gammalib.GModelSpectralPlaw(1.0e-18, -2.0,
                       gammalib.GEnergy(1.0, 'TeV'))
            model    = gammalib.GModelSky(spatial, spectral)
            model.name('Dummy')
            self._models.append(model)
            self._srcname = 'Dummy'
            self['use_model_bkg'].boolean(False)

        # Return
        return

    def _compute_region_separation(self, pnt_dir):
        """
        Compute the separation angle for reflected off regions in radians

        Returns
        -------
        angle : float
            Separation angle of two off regions (radians)
        """
        # Initialise the result
        separation = -1.0

        # Compute offset of reflected regions to pointing position
        offset = pnt_dir.dist_deg(self._src_dir)

        # If shape is a circle then compute apparent diameter of the circle
        # as separation
        if self._srcshape == 'CIRCLE':
            separation = 2.0 * self._rad / offset

        # ... otherwise if shape is a rectangle then compute the opening angle
        # towards combinations of rectangle corners. This method overestimates
        # the real need of space between the ectangles, so the method may be
        # optimised to gain more off regions! Anyway, it is assured that the
        # off regions will never overlap.
        elif self._srcshape == 'RECT':

            # Get the sky directions of the corners of the rectangle
            cs = [self._src_reg[0].corner(icorner) for icorner in range(4)]

            # Compute the 6 opening angles
            combinations = [[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]]
            angles       = [self._compute_posang(pnt_dir, cs[i], cs[j]) \
                            for i,j in combinations]

            # The desired separation is the maximum opening angle
            separation = max(angles) * gammalib.deg2rad

        # Return
        return separation

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

        # Skip observation if it is too close to source
        if self._src_reg.contains(pnt_dir):
            msg = ' Skip because observation is pointed at %.3f deg from source'\
                  % (offset)
            if self._srcshape == 'CIRCLE':
                msg += ' (circle rad=%.3f).' % (self._rad)
            self._log_string(gammalib.NORMAL, msg)
        # ... otherwise
        else:
            posang = pnt_dir.posang_deg(self._src_dir)
            if (self._srcshape == 'CIRCLE') or (self._srcshape == 'RECT'):

                # Determine number of background regions to skip
                N_skip  = self['bkgregskip'].integer()
                N_lim   = 1 + 2 * N_skip

                # Compute the angular separation of reflected regions wrt
                # camera center. The factor 1.05 ensures background regions
                # do not overlap due to numerical precision issues
                alpha = 1.05 * self._compute_region_separation(pnt_dir)

                # Compute number of reflected regions by dividing the angular
                # separation by 2 pi.
                N = int(2.0 * math.pi / alpha)

                # If there are not enough reflected regions then skip the
                # observation ...
                if N < self['bkgregmin'].integer() + N_lim:
                    msg = ' Skip because the number %d of reflected regions '\
                          'for background estimation is smaller than '\
                          '"bkgregmin"=%d.' % (N-N_lim, self['bkgregmin'].integer())
                    self._log_string(gammalib.NORMAL, msg)

                # ... otherwise loop over position angle to create reflected
                # regions
                else:

                    # Append reflected regions
                    alpha    = 360.0 / N
                    dphi_max = 360.0 - alpha * (1 + N_skip)
                    dphi     = alpha         * (1 + N_skip)
                    while dphi <= dphi_max:
                        ctr_dir = pnt_dir.clone()
                        ctr_dir.rotate_deg(posang + dphi, offset)
                        if self._srcshape == 'CIRCLE':
                            region = gammalib.GSkyRegionCircle(ctr_dir, self._rad)
                        elif self._srcshape == 'RECT':
                            # Adjust the posang of the rectangle correspondingly
                            region = gammalib.GSkyRegionRectangle(ctr_dir,
                                                                  self._reg_width,
                                                                  self._reg_height,
                                                                  self._reg_posang +          dphi)
                        if self._has_exclusion:
                            if self._excl_reg.overlaps(region):

                                # Signal region overlap
                                msg = ' Reflected region overlaps with '\
                                      'exclusion region.'
                                self._log_string(gammalib.EXPLICIT, msg)

                                # If region overlaps with exclusion region
                                # try to increment by 10% of angular step
                                dphi += 0.1 * alpha

                            else:
                                regions.append(region)
                                dphi += alpha
                        else:
                            regions.append(region)
                            dphi += alpha

                    # Check again that we have enough background regions
                    # now that we have checked for overlap with exclusion region
                    if regions.size() >= self['bkgregmin'].integer():
                        # Log number of reflected regions
                        msg = ' Use %d reflected regions.' % (regions.size())
                        self._log_string(gammalib.NORMAL, msg)
                    # Otherwise log observation skipped and return empty region container
                    else:
                        msg = ' Skip because the number %d of regions' \
                              'for background estimation not overlapping ' \
                              'with the exclusion region is smaller than ' \
                              '"bkgregmin"=%d.' % \
                              (regions.size(), self['bkgregmin'].integer())
                        self._log_string(gammalib.NORMAL, msg)
                        regions = gammalib.GSkyRegions()

        # Return reflected regions
        return regions

    def _instrument_regions(self, obs, obs_off):
        """
        Compute background regions for Off observation
        
        Calculate background region in Off observation that corresponds to the
        source region in the On observation in instrument coordinates

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation()`
            On CTA observation
        obs_off : `~gammalib.GCTAObservation()`
            Off CTA observation

        Returns
        -------
        regions : `~gammalib.GSkyRegions`
            Container with background region
        """
        # Initialise region container
        regions = gammalib.GSkyRegions()

        # Convert source position in On observation to instrument coordinates
        instdir = obs.pointing().instdir(self._src_dir)

        # Convert instrument position to sky direction for Off observation
        off_dir = obs_off.pointing().skydir(instdir)

        # Build region according to shape specified by user
        # If circle
        if self._srcshape == 'CIRCLE':
            region = gammalib.GSkyRegionCircle(off_dir, self._rad)

        # ... otherwise if rectangle
        elif self._srcshape == 'RECT':
            # Instrument coordinates take sky direction as reference
            # so no need to change the position angle
            region = gammalib.GSkyRegionRectangle(off_dir,
                                                  self._reg_width,
                                                  self._reg_height,
                                                  self._reg_posang)

        # Check if background region overlaps with exclusion region
        is_valid = True
        if self._has_exclusion:
            if self._excl_reg.overlaps(region):
                # Signal region overlap
                msg = ' Background region overlaps with exclusion region.'
                self._log_string(gammalib.EXPLICIT, msg)
                is_valid = False

        # If region is valid then append it to container
        if is_valid:
            regions.append(region)

        # Return
        return regions

    def _set_models(self, results):
        """
        Set models for On/Off fitting

        The method does the following
        - append "OnOff" to the instrument name of all background models
        - fix all spatial and temporal parameters

        Parameters
        ----------
        results : list of dict
            Result dictionaries

        Returns
        -------
        models : `~gammalib.GModels()`
            Model container
        """
        # Write header
        self._log_header1(gammalib.NORMAL, 'Set models')

        # Initialise model container
        models = gammalib.GModels()

        # Initialies stacked model flag
        has_stacked_model = False

        # Loop over all models in observation and append "OnOff" to instrument
        # names
        for model in self._models:

            # Initialise model usage
            use_model = False

            # If model is a background model then check if it will be
            # used
            if 'GCTA' in model.classname():

                # Skip model if background model should not be used
                if not self['use_model_bkg'].boolean():
                    self._log_string(gammalib.NORMAL, ' Skip "%s" model "%s" (%s)' % \
                         (model.instruments(), model.name(), model.ids()))
                    continue

                # Check if model corresponds to one of the relevant
                # observations
                for result in results:
                    if model.is_valid(result['instrument'], result['id']):
                        if result['bkg_reg'].size() > 0:
                            use_model = True
                        break

                # If stacked analysis is requested then just use for model
                # and remove instrument ID
                if self['stack'].boolean():

                    # If there is already a model for stacked analysis then
                    # skip this one
                    if has_stacked_model:
                        msg = ' Skip "%s" model "%s" (%s). There is already ' \
                              'a model for stacked analysis.' % \
                              (model.instruments(), model.name(), model.ids())
                        self._log_string(gammalib.NORMAL, msg)
                        continue

                    # ... otherwise use model for stacked analysis
                    else:
                        has_stacked_model = True
                        use_model         = True
                        model.ids('')

                # Append "OnOff" to instrument name
                model.instruments(model.instruments()+'OnOff')

            # ... otherwise, if model is not a background model then use it
            else:
                use_model = True

            # If model is relevant then append it now to the model
            # container
            if use_model:

                # Log model usage
                self._log_string(gammalib.NORMAL, ' Use "%s" model "%s" (%s)' % \
                     (model.instruments(), model.name(), model.ids()))

                # Append model to container
                models.append(model)

            # ... otherwise signal that model is skipped
            else:
                self._log_string(gammalib.NORMAL, ' Skip "%s" model "%s" (%s)' % \
                     (model.instruments(), model.name(), model.ids()))

        # Return model container
        return models

    def _set_statistic(self, obs):
        """
        Set statistic for observation

        If the "use_model_bkg" is true then set statistic to "cstat",
        otherwise set it to "wstat"

        Parameters
        ----------
        obs : `~gammalib.GObservation()`
            Observation

        Returns
        -------
        obs : `~gammalib.GObservation()`
            Observation
        """
        # Set statistic according to background model usage
        if self['use_model_bkg'].boolean():
            obs.statistic('cstat')
        else:
            obs.statistic('wstat')

        # Return observation
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

    def _set_background_regions(self, obs, obs_off=None):
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

        # If reflected background is requested then create reflected
        # background regions
        if self['bkgmethod'].string() == 'REFLECTED':
            bkg_reg = self._reflected_regions(obs)

        # ... otherwise if custom background is requested then get the
        # background regions from the observation. We use a copy here since
        # otherwise the background regions go out of scope once the observations
        # are replaced by the On/Off observations.
        elif self['bkgmethod'].string() == 'CUSTOM':
            bkg_reg = obs.off_regions().copy()

        # ... otherwise if dedicated Off runs are use then
        # use background region that correspond to the same instrument coordinates
        if self['bkgmethod'].string() == 'OFF':
            bkg_reg = self._instrument_regions(obs,obs_off)

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

        # Retrieve dedicated Off observation if it exists
        if not self.obs_off().is_empty():
            obs_off = self.obs_off()[i]
        # Otherwise use the same as On
        else:
            obs_off = self.obs()[i]

        # Log header
        self._log_header3(gammalib.NORMAL,'%s observation "%s"' % \
                          (obs.instrument(), obs.id()))

        # Skip non CTA observations
        if obs.classname() != 'GCTAObservation':
            self._log_string(gammalib.NORMAL, ' Skip because not a "GCTAObservation"')

        # Otherwise calculate On/Off spectra
        else:
            # Get background model usage flag and log flag
            use_model_bkg = self['use_model_bkg'].boolean()
            if use_model_bkg:
                msg = ' Use background model.'
            else:
                msg = ' Background model not used, assume constant backround rate.'
            self._log_string(gammalib.NORMAL, msg)

            # Get offset angle of source
            pnt_dir = obs.pointing().dir()
            offset = pnt_dir.dist_deg(self._src_dir)

            # Skip observation if it is pointed too far from the source
            if offset >= self['maxoffset'].real():
                msg = ' Skip because observation is pointed at %.3f deg >= ' \
                      '"maxoffset=%.3f" from source.' \
                      % (offset, self['maxoffset'].real())
                self._log_string(gammalib.NORMAL, msg)
            # ... otherwise continue to process
            else:

                # Set background regions for this observation
                bkg_reg = self._set_background_regions(obs,obs_off)

                # If there are any background regions then create On/Off observation
                # and append it to the output container
                if bkg_reg.size() >= 0:

                    # Create On/Off observation
                    onoff = gammalib.GCTAOnOffObservation(obs, obs_off,
                                                          self._models,
                                                          self._srcname,
                                                          self._etruebounds,
                                                          self._ebounds,
                                                          self._src_reg,
                                                          bkg_reg,
                                                          use_model_bkg)

                    # Set On/Off observation ID
                    onoff.id(obs.id())

                # Otherwise log observation skipped
                else:
                    msg = ' Skip because no valid Off regions could be determined'
                    self._log_string(gammalib.NORMAL, msg)

        # Construct dictionary with results
        result = {'onoff'     : onoff,
                  'bkg_reg'   : bkg_reg,
                  'instrument': obs.instrument(),
                  'id'        : obs.id()}

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

                # Set statistic according to background model usage
                obs = self._set_statistic(result['onoff'])

                # Append observation to observation container
                outobs.append(obs)

                # Append background regions
                self._bkg_regs.append({'regions': result['bkg_reg'],
                                      'id': result['id']})

        # Return observation container
        return outobs


    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Observation')
        if not self.obs_off().is_empty():
            self._log_observations(gammalib.NORMAL, self._obs_off, 'Off Observation')

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
        results        = []

        # If there is more than one observation and we use multiprocessing
        if self._nthreads > 1 and self.obs().size() > 1:

            # Compute observations
            args        = [(self, '_process_observation', i)
                           for i in range(self.obs().size())]
            poolresults = mputils.process(self._nthreads, mputils.mpfunc, args)

            # Construct results
            for i in range(self.obs().size()):
                result = poolresults[i][0]
                outobs = self._unpack_result(outobs, result)
                results.append(result)
                self._log_string(gammalib.TERSE, poolresults[i][1]['log'], False)

        # Otherwise, loop through observations and generate pha, arf, rmf files
        else:
            for i in range(self.obs().size()):
                # Process individual observation
                result = self._process_observation(i)
                outobs = self._unpack_result(outobs, result)
                results.append(result)

        # Stack observations
        if outobs.size() > 1 and self['stack'].boolean():

            # Write header
            self._log_header1(gammalib.NORMAL, 'Stacking %d observations' %
                              (outobs.size()))

            # Stack observations
            stacked_obs = gammalib.GCTAOnOffObservation(outobs)

            # Set statistic according to background model usage
            stacked_obs = self._set_statistic(stacked_obs)

            # Put stacked observations in output container
            outobs = gammalib.GObservations()
            outobs.append(stacked_obs)

        # Create models that allow On/Off fitting
        models = self._set_models(results)

        # Set models in output container
        outobs.models(models)

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
        outobs   = self['outobs'].filename()
        outmodel = self['outmodel'].filename()
        prefix   = self['prefix'].string()
        clobber  = self['clobber'].boolean()

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

            # Stamp files
            self._stamp(onname)
            self._stamp(offname)
            self._stamp(arfname)
            self._stamp(rmfname)

            # Log file names
            self._log_value(gammalib.NORMAL, 'PHA on file', onname)
            self._log_value(gammalib.NORMAL, 'PHA off file', offname)
            self._log_value(gammalib.NORMAL, 'ARF file', arfname)
            self._log_value(gammalib.NORMAL, 'RMF file', rmfname)

        # Save observation definition XML file
        self.obs().save(outobs)

        # Save model definition XML file
        self.obs().models().save(outmodel)

        # Log file names
        self._log_value(gammalib.NORMAL, 'Obs. definition XML file', outobs.url())
        self._log_value(gammalib.NORMAL, 'Model definition XML file', outmodel.url())

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

    def exclusion_map(self, object=None):
        """
        Return and optionally set the exclusion regions map

        Parameters
        ----------
        object : `~gammalib.GSkyRegion` or `~gammalib.GSkyMap` or `~gammalib.GFilename`
            Exclusion regions object

        Returns
        -------
        region : `~gammalib.GSkyRegionMap`
            Exclusion regions map
        """
        # If a regions object is provided then set the regions ...
        if object is not None:
            self._excl_reg = gammalib.GSkyRegionMap(object)

        # Return
        return self._excl_reg

    def obs_off(self, obs=None):
        """
        Return and optionally set the Off observations

        Parameters
        ----------
        obs : `~gammalib.GCTAObservations`
            Off observations container

        Returns
        -------
        observation container : `~gammalib.GCTAObservations`
            Off observations container
        """
        # If an observation container is provided then set the Off observations ...
        if obs is not None:
            self._obs_off = obs

        # Return
        return self._obs_off


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csphagen(sys.argv)

    # Execute application
    app.execute()
