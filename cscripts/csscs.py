#! /usr/bin/env python
# ==========================================================================
# Spectral component separation script
#
# Copyright (C) 2020-2022 Luigi Tibaldo
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
from cscripts import mputils
from cscripts import obsutils


# =========== #
# csscs class #
# =========== #
class csscs(ctools.csobservation):
    """
    Spectral component separation script
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor

        Parameters
        ----------
        argv : list of str
            List of IRAF command line parameter strings of the form
            ``parameter=3``.
        """
        # Initialise application by calling the base class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise data members
        self._nthreads      = 0
        self._srcnames      = []
        self._method        = None
        self._fits          = None
        self._excl_reg_map  = None  # Exclusion region map for on/off analysis
        self._excl_reg_name = None

        # Return
        return

    # State methods for pickling
    def __getstate__(self):
        """
        Extend ctools.csobservation __getstate__ method

        Returns
        -------
        state : dict
            Pickled instance
        """
        # Set pickled dictionary
        state = {'base':          ctools.csobservation.__getstate__(self),
                 'nthreads':      self._nthreads,
                 'srcnames':      self._srcnames,
                 'method':        self._method,
                 'fits':          self._fits,
                 'excl_reg_map':  self._excl_reg_map,
                 'excl_reg_name': self._excl_reg_name}

        # Return pickled dictionary
        return state

    def __setstate__(self, state):
        """
        Extend ctools.csobservation __setstate__ method

        Parameters
        ----------
        state : dict
            Pickled instance
        """
        # Set state
        ctools.csobservation.__setstate__(self, state['base'])
        self._nthreads      = state['nthreads']
        self._srcnames      = state['srcnames']
        self._method        = state['method']
        self._fits          = state['fits']
        self._excl_reg_map  = state['excl_reg_map']
        self._excl_reg_name = state['excl_reg_name']

        # Return
        return

    # Private methods
    def _get_onoff_parameters(self):
        """
        Get On/Off analysis parameters.

        This is done here rather than in ctools.is_onoff because the exclusion
        map needs to be handled differently as an automatic parameter but
        queried only if not already set and we need to verify is not empty.
        Also many parameters are already queried for all methods

        TODO: verify if this can be made more uniform with other scripts
        """
        # Exclusion map
        if (self._excl_reg_map is not None) and (self._excl_reg_map.map().npix() > 0):
            # Exclusion map set and is not empty
            pass
        elif self['inexclusion'].is_valid():
            self._excl_reg_name = self['inexclusion'].filename()
            self._excl_reg_map = gammalib.GSkyRegionMap(self._excl_reg_name)
        else:
            msg = 'csscs in On/Off mode requires input exclusion region.'
            raise RuntimeError(msg)

        # Other csphagen parameters
        self['enumbins'].integer()
        if self['bkgmethod'].string() == 'REFLECTED':
            self['bkgregmin'].integer()
            self['bkgregskip'].integer()
        self['maxoffset'].real()
        self['etruemin'].real()
        self['etruemax'].real()
        self['etruebins'].integer()

        # Return
        return

    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Set observation if not done before
        if self.obs().is_empty():
            self.obs(self._get_observations())

        # Set observation statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

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

        # Check if there are On/Off or non-IACT observations
        if n_onoff > 0 or n_other > 0:
            msg = 'On/Off or non-CTA observations found. ' \
                  'csscs does not support this type of observations.'
            raise RuntimeError(msg)

        # ... otherwise check if there is a mix of binned and unbinned
        elif n_unbinned > 0 and n_binned > 0:
            msg = 'Mix of unbinned and binned CTA observations ' \
                  'found in observation container. csscs does not ' \
                  'support this mix.'
            raise RuntimeError(msg)

        # Set models if there are none in the container
        if self.obs().models().size() == 0:
            self.obs().models(self['inmodel'].filename())

        # Query the map definition parameters
        self['xref'].real()
        self['yref'].real()
        self['coordsys'].string()
        self['proj'].string()
        self['nxpix'].integer()
        self['nypix'].integer()
        self['binsz'].real()

        # Query radius of ROI for component separation
        self['rad'].real()

        # Query energy boundaries
        self['emin'].real()
        self['emax'].real()

        # If we have unbinned observations query analysis method
        if n_unbinned > 0:
            self._method = self['method'].string()
            # If method is Onoff query csphagen parameters
            if self._method == 'ONOFF':
                self._get_onoff_parameters()
        # Otherwise set method to binned
        else:
            self._method = 'BINNED'

        # Query target source names
        srcnames = self['srcnames'].string()

        # Fashion target sources into Python list. Strip leading and trailing
        # spaces from names.
        self._srcnames = srcnames.split(';')
        for s, name in enumerate(self._srcnames):
            self._srcnames[s] = name.strip()

        # Verify that the target sources are in the model
        src_in_model = True
        for name in self._srcnames:
            src_in_model = src_in_model and self.obs().models().contains(name)
        if src_in_model == False:
            msg = 'Not all target sources are present in input model.'
            raise RuntimeError(msg)

        # For On/Off analysis check the number of gamma-ray sources in the model
        if self._method == 'ONOFF':

            # Initialise number of sources
            nsources = 0

            # Loop over models and verify if model is of type GModelSky
            for model in self.obs().models():
                if model.classname() == 'GModelSky':
                    nsources += 1

            # If there are background gamma-ray sources in the model
            # throw runtime error
            if nsources > len(self._srcnames):
                msg = 'Background gamma-ray sources found in the model. ' \
                      'On/Off analysis does not support this feature.'
                raise RuntimeError(msg)

            # ... otherwise continue
            else:
                pass

        # Query other hidden parameters, just in case
        self['edisp'].boolean()
        self['calc_ulim'].boolean()
        self['calc_ts'].boolean()
        self['fix_bkg'].boolean()
        self['fix_srcs'].boolean()

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].query()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Set number of processes for multiprocessing
        self._nthreads = mputils.nthreads(self)

        # Return
        return

    def _create_map(self):
        """
        Create a gammalib.GSkyMap object to store the output

        Returns
        -------
        skymap : `~gammalib.GSkyMap`
            Sky map
        """
        # Create sky map
        skymap = gammalib.GSkyMap(self['proj'].string(), self['coordsys'].string(),
                                  self['xref'].real(), self['yref'].real(),
                                  -self['binsz'].real(), self['binsz'].real(),
                                  self['nxpix'].integer(), self['nypix'].integer(),
                                  1)

        # Return sky map
        return skymap

    def _adjust_model_pars(self):
        """
        Adjust model parameters
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Adjust model parameters')

        # Adjust model parameters based on input user parameters
        for model in self.obs().models():

            # Initialise TS flag for all models to false
            model.tscalc(False)

            # Log model name
            self._log_header3(gammalib.EXPLICIT, model.name())

            # In On/Off mode we cannot support multiple source morphologies
            # Thus, if we have more than one target we will set them all
            # to isotropic
            if self._method == 'ONOFF' and len(self._srcnames) > 1:

                # Check if it is a source
                try:
                    gammalib.GModelSky(model)
                    # In this case check if the spatial model is already isotropic
                    if model.spatial() == 'DiffuseIsotropic':
                        pass

                    # ... otherwise change it to isotropic
                    else:
                        msg = ' Setting spatial model to diffuse isotropic'
                        self._log_string(gammalib.EXPLICIT, msg)
                        model.spatial(gammalib.GModelSpatialDiffuseConst())
                except:
                    pass

            # Freeze all parameters except the normalization
            # which is assumed to be the first spectral parameter
            normpar = model.spectral()[0]
            for par in model:
                if par.name() == normpar.name():
                    pass
                else:
                    if par.is_free():
                        self._log_string(gammalib.EXPLICIT, ' Fixing "' + par.name() + '"')
                        par.fix()

            # Deal with the target sources
            if model.name() in self._srcnames:

                # Free the normalisation parameter which is assumed to be
                # the first spectral parameter
                if normpar.is_fixed():
                    self._log_string(gammalib.EXPLICIT, ' Freeing "' + normpar.name() + '"')
                    normpar.free()

                # Optionally compute Test Statistic value
                if self['calc_ts'].boolean():
                    model.tscalc(True)

            # Deal with background models
            elif self['fix_bkg'].boolean() and not model.classname() == 'GModelSky':

                if normpar.is_free():
                    self._log_string(gammalib.EXPLICIT, ' Fixing "' + normpar.name() + '"')
                    normpar.fix()


            # Deal with background sources
            elif self['fix_srcs'].boolean() and model.classname() == 'GModelSky':

                if normpar.is_free():
                    self._log_string(gammalib.EXPLICIT, ' Fixing "' + normpar.name() + '"')
                    normpar.fix()

        # Return
        return

    def _mask_cube(self, obs, ra, dec, rad):
        """
        Mask cube observation

        Parameters
        ----------
        obs : `gammalib.GCTAObservation`
            Input observation of type cube
        ra : float
            Right Ascension (deg)
        dec : float
            Declination (deg)
        rad : float
            Radius (deg)

        Returns
        -------
        new_obs : `gammalib.GCTAObservations`
            Output observations with masked cubes
        """
        # Filter cube according to ROI and user energy range
        cubemask = ctools.ctcubemask(obs)
        cubemask['ra'] = ra
        cubemask['dec'] = dec
        cubemask['rad'] = rad
        cubemask['regfile'] = 'NONE'
        cubemask['emin'] = self['emin'].real()
        cubemask['emax'] = self['emax'].real()

        # If chatter level is verbose and debugging is requested then
        # switch also on the debug model in ctcubemask
        if self._logVerbose() and self._logDebug():
            cubemask['debug'] = True

        # Mask cube
        cubemask.run()

        # Get deep copy of filtered observations
        new_obs = cubemask.obs().copy()

        # Return
        return new_obs

    def _mask_evlist(self, obs, ra, dec, rad):
        """
        Mask event list

        Parameters
        ----------
        obs : `gammalib.GCTAObservations`
            Input observations of type event list
        ra : float
            Right Ascension (deg)
        dec : float
            Declination (deg)
        rad : float
            Radius (deg)

        Returns
        -------
        new_obs : `gammalib.GCTAObservations`
            Output observations with masked event lists
        """
        # Filter event list according to ROI and user energy range
        select = ctools.ctselect(obs)
        select['usepnt'] = False
        select['ra']     = ra
        select['dec']    = dec
        select['rad']    = rad
        select['emin']   = self['emin'].real()
        select['emax']   = self['emax'].real()
        select['tmin']   = 'INDEF'
        select['tmax']   = 'INDEF'

        # If chatter level is verbose and debugging is requested then
        # switch also on the debug model in ctcubemask
        if self._logVerbose() and self._logDebug():
            select['debug'] = True

        # Select event list
        select.run()

        # Get deep copy of filtered observation
        new_obs = select.obs().copy()

        # Return
        return new_obs

    def _mask_observations(self, ra, dec, rad):
        """
        Create observations restricted to circular ROI

        Parameters
        ----------
        ra : float
            Right Ascension (deg)
        dec : float
            Declination (deg)
        rad : float
            Radius (deg)

        Returns
        -------
        new_obs : `~gammalib.GObservations`
            Observations in circular ROI
        """
        # Determine type of masking according to observation type and analysis
        # method
        if self._method == 'UNBINNED':
            new_obs = self._mask_evlist(self.obs(), ra, dec, rad)
        elif self._method == 'BINNED':
            new_obs = self._mask_cube(self.obs(), ra, dec, rad)
        elif self._method == 'ONOFF':
            new_obs = obsutils.get_onoff_obs(self, self.obs(), nthreads=1,
                                             ra=ra, dec=dec,
                                             srcname=self._srcnames[0])

        # Return
        return new_obs

    def _pixel_analysis(self, inx):
        """
        Performs analysis over the region of interest corresponding to a
        single pixel of output map

        Parameters
        ----------
        inx : int
            Pixel index

        Returns
        -------
        result : dict
            Results
        """

        # Create template for output map to extract pixel coordinates
        outmap = self._create_map()

        # Determine pixel centre coordinates
        pixdir = outmap.inx2dir(inx)
        ra     = pixdir.ra_deg()
        dec    = pixdir.dec_deg()

        # Write header for spatial bin
        msg = 'Spatial bin (%d) centred on (R.A.,Dec) = (%f,%f) deg' % (inx, ra, dec)
        self._log_header2(gammalib.TERSE, msg)

        # Initialize results
        result = {}
        for name in self._srcnames:
            result[name] = {'flux': 0.0,
                            'flux_err': 0.0,
                            'TS': 0.0,
                            'ulimit': -1.0}

        # Mask observations
        self._log_header3(gammalib.NORMAL, 'Masking observations')
        masked_obs = self._mask_observations(ra, dec, self['rad'].real())

        # Continue only if we have at least one observation
        if masked_obs.size() > 0:

            # Set up likelihood analysis
            self._log_header3(gammalib.NORMAL, 'Performing fit in region')
            like = ctools.ctlike(masked_obs)
            like['edisp']    = self['edisp'].boolean()
            like['nthreads'] = 1  # Avoids OpenMP conflict

            # If chatter level is verbose and debugging is requested then
            # switch also on the debug model in ctlike
            if self._logVerbose() and self._logDebug():
                like['debug'] = True

            # Perform maximum likelihood fit
            like.run()

            # Write model results for explicit chatter level
            self._log_string(gammalib.EXPLICIT, str(like.obs().models()))

            # Continue only if log-likelihood is non-zero
            if like.obs().logL() != 0.0:

                # Prepare objects for flux extraction
                # ROI
                centre = gammalib.GSkyDir()
                centre.radec_deg(ra, dec)
                roi = gammalib.GSkyRegionCircle(centre, self['rad'].real())

                # Energy boundaries
                emin = gammalib.GEnergy(self['emin'].real(), 'TeV')
                emax = gammalib.GEnergy(self['emax'].real(), 'TeV')

                # Loop over target sources
                for name in self._srcnames:

                    # Get source
                    source = like.obs().models()[name]

                    # Get flux
                    # Integrate spectral model between emin and emax
                    flux = source.spectral().flux(emin, emax)

                    # Calculate correction factor
                    # Spatial model flux over ROI divided by ROI solid angle
                    corr_factor  = source.spatial().flux(roi)
                    corr_factor /= gammalib.twopi * (1.0 - math.cos(math.radians(self['rad'].real())))

                    # Multiply flux by correction factor
                    flux                *= corr_factor
                    result[name]['flux'] = flux

                    # Get flux error
                    # Normalization parameter
                    normpar = source.spectral()[0]

                    # Only normalization parameter free
                    # Relative error on flux is same as on normalization
                    # Avoid zero division error
                    if normpar.value() > 0.:
                        flux_error = flux * normpar.error() / normpar.value()
                    else:
                        flux_error = 0.
                    result[name]['flux_error'] = flux_error

                    # If requested get TS
                    if self['calc_ts'].boolean():
                        result[name]['TS'] = source.ts()

                    # If requested compute upper flux limit
                    if self['calc_ulim'].boolean():

                        # Logging information
                        self._log_header3(gammalib.NORMAL,
                                          'Computing upper limit for source ' + name)

                        # Create upper limit object
                        ulimit = ctools.ctulimit(like.obs())
                        ulimit['srcname'] = name
                        ulimit['emin'] = self['emin'].real()
                        ulimit['emax'] = self['emax'].real()

                        # If chatter level is verbose and debugging is requested
                        # then switch also on the debug model in ctulimit
                        if self._logVerbose() and self._logDebug():
                            ulimit['debug'] = True

                        # Try to run upper limit and catch exceptions
                        try:
                            ulimit.run()
                            ulimit_value = ulimit.flux_ulimit()
                            # Multiply by correction factor to get flux per solid angle in ROI
                            ulimit_value *= corr_factor
                            result[name]['ulimit'] = ulimit_value
                        except:
                            self._log_string(gammalib.NORMAL, 'Upper limit '
                                                              'calculation failed.')

            # if likelihood is zero
            else:
                value = 'Likelihood is zero. Bin is skipped.'
                self._log_value(gammalib.TERSE, '(R.A.,Dec) = (%f,%f) deg' % (ra, dec), value)

        # if observation size = 0
        else:
            value = 'Size of masked observations is zero. Bin is skipped.'
            self._log_value(gammalib.TERSE, '(R.A.,Dec) = (%f,%f) deg' % (ra, dec), value)

        # Return result
        return result

    def _write_hdu_keywords(self, hdu, is_flux=False):
        """
        Append cards to HDU

        Parameters
        ----------
        hdu : `~gammalib.GFitsHDU`
            HDU
        """
        # Flux units
        if is_flux:
            hdu.card('BUNIT', 'ph/cm2/s/sr', 'Photon flux')

        # Minimum and maximum energy
        hdu.card('E_MIN', self['emin'].real(), '[TeV] Lower energy boundary')
        hdu.card('E_MAX', self['emax'].real(), '[TeV] Upper energy boundary')
        hdu.card('EUNIT', 'TeV', 'Units for E_MIN and E_MAX')

        # ROI size
        hdu.card('ROISZ', self['rad'].real(),
                 '[deg] Size of ROI used for component separation')

        # Analysis method
        hdu.card('METHOD', self._method, 'Analysis method')

        # For On/Off log information
        if self._method == 'ONOFF':

            # Exclusion map
            if self._excl_reg_name is not None:
                hdu.card('EXCLMAP', self._excl_reg_name.url(), 'Exclusion map name')

            # Model background?
            hdu.card('BKGMODEL', self['use_model_bkg'].boolean(), 'Use background model?')

            # Fit statistic
            hdu.card('STAT', self['statistic'].string(), 'Fit statistic')

        # Return
        return

    def _fill_fits(self, results):
        """
        Fill FITS object to store the results

        Parameters
        ----------
        result : `dict`
            Results

        Returns
        -------
        fits : `~gammalib.GFits`
            FITS object
        """
        # Create new Fits object
        fits = gammalib.GFits()

        # Append empty primary HDU
        fits.append(gammalib.GFitsImageDouble())

        # Loop over target sources
        for s, name in enumerate(self._srcnames):

            # Create skymaps to store results
            # Maps for TS and upper limits created and filled to avoid if statements
            # Will not be saved in the end if not requested
            fmap   = self._create_map()
            errmap = self._create_map()
            tsmap  = self._create_map()
            ulmap  = self._create_map()

            # Loop over pixels and fill maps
            for inx in range(fmap.npix()):
                fmap[inx, 0]   = results[inx][name]['flux']
                errmap[inx, 0] = results[inx][name]['flux_error']
                tsmap[inx, 0]  = results[inx][name]['TS']
                ulmap[inx, 0]  = results[inx][name]['ulimit']

            # Uppercase name to adhere to gammalib convention
            Name = gammalib.toupper(name)

            # Write maps to Fits
            fhdu = fmap.write(fits, Name + ' FLUX')
            self._write_hdu_keywords(fhdu, is_flux=True)
            errhdu = errmap.write(fits, Name + ' FLUX ERROR')
            self._write_hdu_keywords(errhdu, is_flux=True)

            # If requested write TS map to fits
            if self['calc_ts'].boolean():
                tshdu = tsmap.write(fits, Name + ' TS')
                self._write_hdu_keywords(tshdu)

            # If requested write upper limit map to fits
            if self['calc_ulim'].boolean():
                ulhdu = ulmap.write(fits, Name + ' FLUX UPPER LIMIT')
                self._write_hdu_keywords(ulhdu, is_flux=True)

        # Return
        return fits

    def _get_skymap(self, name, quantity):
        """
        Fetches skymap from FITS container

        Parameters
        ----------
        name : str
            Source name
        quantity : str
            Quantity

        Returns
        -------
        skymap : `~gammalib.GSkyMap`
            Skymap

        Raises
        ------
        NameError
        """

        # Initialise skymap object
        skymap = None

        # Proceed only if fits has been filled
        if self._fits != None:

            # Check that the name is in the list of target sources
            if name in self._srcnames:

                # Uppercase name
                Name = gammalib.toupper(name)

                # Assemble HDU name
                hduname = Name + ' ' + quantity

                # Try to fetch HDU
                try:
                    hdu = self._fits[hduname]
                    # Convert to GSkyMap
                    skymap = gammalib.GSkyMap(hdu)
                # Catch exception if HDU does not exist
                except:
                    msg = 'HDU "' + hduname + '" not found in FITS container.'
                    raise RuntimeError(msg)

            # ... otherwise throw name error
            else:
                print('ERROR: Source "' + name + '" not found in list of target sources.')
                raise NameError(name)

        # Return
        return skymap

    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Write input observation container into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Input observation')

        # Adjust model parameters based on input user parameters
        self._adjust_model_pars()

        # Compute number of pixels of output map
        npix = self._create_map().npix()

        # If more than a single thread is requested then use multiprocessing
        if self._nthreads > 1:

            # Compute values in pixels
            args = [(self, '_pixel_analysis', i)
                    for i in range(npix)]
            poolresults = mputils.process(self._nthreads, mputils.mpfunc, args)

            # Construct results
            results = []
            for i in range(npix):
                results.append(poolresults[i][0])
                self._log_string(gammalib.TERSE, poolresults[i][1]['log'], False)

        # Otherwise loop over pixels and run the pixel analysis
        else:
            results = []
            for i in range(npix):
                results.append(self._pixel_analysis(i))

        # Fill output Fits
        self._fits = self._fill_fits(results)

        # Stamp FITS file
        self._stamp(self._fits)

        # Return
        return

    def save(self):
        """
        Save maps to Fits flie
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save maps')

        # Continue only if FITS file is valid
        if self._fits != None:
            # Get outmap parameter
            outfile = self['outfile'].filename()

            # Log file name
            self._log_value(gammalib.NORMAL, 'Output file', outfile.url())

            # Save spectrum
            self._fits.saveto(outfile, self['clobber'].boolean())

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
        # If a regions object is provided then set the regions
        if object is not None:
            self._excl_reg_map = gammalib.GSkyRegionMap(object)

        # Return
        return self._excl_reg_map

    def fits(self):
        """
        Return copy of the fits container

        Returns
        -------
        fits : `~gammalib.GFits`
            FITS file containing the maps
        """
        # Return
        return self._fits.copy()

    def flux(self, name):
        """
        Return flux skymap

        Parameters
        ----------
        name : str
            Source name

        Returns
        -------
        skymap : `~gammalib.GSkyMap`
            Flux skymap
        """
        # Get flux skymap
        skymap = self._get_skymap(name, 'FLUX')

        # Return flux sky map
        return skymap

    def flux_error(self, name):
        """
        Return flux error skymap

        Parameters
        ----------
        name : str
            Source name

        Returns
        -------
        skymap : `~gammalib.GSkyMap`
            Flux error skymap
        """
        # Get flux error sky map
        skymap = self._get_skymap(name, 'FLUX ERROR')

        # Return flux error sky map
        return skymap

    def ts(self, name):
        """
        Return TS skymap

        Parameters
        ----------
        name : str
            Source name

        Returns
        -------
        skymap : `~gammalib.GSkyMap`
            TS skymap
        """
        # Check that TS computation is requested
        if self['calc_ts'].boolean():
            skymap = self._get_skymap(name, 'TS')

        # ... otherwise raise error
        else:
            msg = 'TS computation not requested. Change user parameter ' \
                  '"calc_ts" to True to obtain TS.'
            raise RuntimeError(msg)

        # Return TS map
        return skymap

    def ulimit(self, name):
        """
        Return flux upper limit skymap

        Parameters
        ----------
        name : str
            Source name

        Returns
        -------
        skymap : `~gammalib.GSkyMap`
            Flux upper limit skymap
        """
        # Check that upper limit computation is requested
        if self['calc_ulim'].boolean():
            skymap = self._get_skymap(name, 'FLUX UPPER LIMIT')

        # ... otherwise raise error
        else:
            msg = 'Upper limit computation not requested. Change user parameter ' \
                  '"calc_ulimit" to True to obtain upper limit.'
            raise RuntimeError(msg)

        # Return upper limit map
        return skymap


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csscs(sys.argv)

    # Execute application
    app.execute()
