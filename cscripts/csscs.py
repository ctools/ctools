#! /usr/bin/env python
# ==========================================================================
# Spectral component separation script
#
# Copyright (C) 2020 Luigi Tibaldo
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
from cscripts import mputils
from cscripts import obsutils
import math


# ================= #
# csscs class #
# ================= #
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
        self._nthreads = 0
        self._roisz = 0.
        self._srcnames = []
        self._fits = None
        self._excl_reg_map = None # Exclusion region map for on/off analysis

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
        # Set pickled dictionary
        state = {'base'     : ctools.csobservation.__getstate__(self),
                 'nthreads' : self._nthreads,
                 'roisz'    : self._roisz,
                 'srcnames' : self._srcnames,
                 'fits'     : self._fits,
                 'excl_reg_map' : self._excl_reg_map}

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
        self._nthreads = state['nthreads']
        self._roisz = state['roisz']
        self._srcnames = state['srcnames']
        self._fits = state['fits']
        self._excl_reg_map = state['excl_reg_map']

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """

        # Set observation if not done before
        if self.obs().is_empty():
            self.obs(self._get_observations())

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

        # Check that we have at least one unbinned or binned CTA observation
        if n_unbinned == 0 and n_binned == 0:
            msg = 'No unbinned or binned CTA observations found. ' \
                  'Please provide at least one unbinned or binned ' \
                  'CTA observation.'
            raise RuntimeError(msg)

        # Set observation statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

        # Set models if there are none in the container
        if self.obs().models().size() == 0:
            self.obs().models(self['inmodel'].filename())

        # Query exclusion region
        self["inexclusion"].query()

        # Query the map definition parameters
        self['xref'].real()
        self['yref'].real()
        self['coordsys'].string()
        self['proj'].string()
        self['nxpix'].integer()
        self['nypix'].integer()

        # Compute the size of the ROIs for the analysis in deg
        self._roisz = self['binsz'].real() * self['roisz'].real()

        # Query energy boundaries
        self["emin"].real()
        self["emax"].real()

        # If we have unbinned observations query the analysis method
        if n_unbinned > 0:
            # If method is ONOFF query some csphagen parameters
            if self['method'].string() == 'ONOFF':
                pass
                # TODO: implement Onoff analysis
                # self["enumbins"].integer()
                # self['use_model_bkg'].boolean()
                # if self["bkgmethod"].string() == "REFLECTED":
                #     self["srcshape"].string()
                #     self["bkgregmin"].integer()
                #     self["maxoffset"].real()
                #     self["etruemin"].real()
                #     self["etruemax"].real()
                #     self["etruebins"].integer()

        # Query target source names
        srcnames = self['srcnames'].string()

        # Fashion target sources into Python list
        self._srcnames = srcnames.split(';')
        # Strip leading and trailing spaces from names
        for s, name in enumerate(self._srcnames):
            self._srcnames[s] = name.strip()

        # Verify that the target sources are in the model
        src_in_model = True
        for name in self._srcnames:
            src_in_model = src_in_model and self.obs().models().contains(name)
        if src_in_model == False:
            msg = 'Not all target sources are present in input model.'
            raise RuntimeError(msg)

        # Query the hidden parameters, just in case
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
        """

        skymap = gammalib.GSkyMap(self['proj'].string(), self['coordsys'].string(),
                                  self['xref'].real(), self['yref'].real(),
                                  -self['binsz'].real(), self['binsz'].real(),
                                  self['nxpix'].integer(), self['nypix'].integer(),
                                  1)

        return skymap

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
                    self._log_string(gammalib.EXPLICIT, ' Freeing "'+normpar.name()+'"')
                    normpar.free()

                # Optionally compute Test Statistic value
                if self['calc_ts'].boolean():
                    model.tscalc(True)

            # Deal with background models
            elif self['fix_bkg'].boolean() and not model.classname() == 'GModelSky':

                if normpar.is_free():
                    self._log_string(gammalib.EXPLICIT, ' Fixing "'+normpar.name()+'"')
                    normpar.fix()


            # Deal with background sources
            elif self['fix_srcs'].boolean() and model.classname() == 'GModelSky':

                if normpar.is_free():
                    self._log_string(gammalib.EXPLICIT, ' Fixing "'+normpar.name()+'"')
                    normpar.fix()

        # Return
        return

    def _mask_cube(self,obs,ra,dec,rad):
        """
        Mask cube observation

        Parameters
        ----------
        obs     : `gammalib.GCTAObservation` Input observation of type cube
        ra      : `float` R.A. (deg)
        dec     : `float` Dec (deg)
        rad     : `float` radius (deg)

        Returns
        -------
        new_obs : `gammalib.GCTAObservation`
            Output observation with masked cube
        """

        # Put observation in container
        container = gammalib.GObservations()
        container.append(obs)

        # Filter cube according to ROI and user energy range
        cubemask = ctools.ctcubemask(container)
        cubemask['ra'] = ra
        cubemask['dec'] = dec
        cubemask['rad'] = rad
        cubemask['regfile'] = 'NONE'
        cubemask['emin'] = self['emin'].real()
        cubemask['emax'] = self['emax'].real()
        cubemask.run()

        # Get deep copy of filtered observation
        new_obs = cubemask.obs()[0].copy()

        # Return
        return new_obs

    def _mask_evlist(self,obs,ra,dec,rad):
        """
        Mask event list

        Parameters
        ----------
        obs     : `gammalib.GCTAObservation` Input observation of type cube
        ra      : `float` R.A. (deg)
        dec     : `float` Dec (deg)
        rad     : `float` radius (deg)

        Returns
        -------
        new_obs : `gammalib.GCTAObservation`
            Output observation with masked event list
        """

        # Put observation in container
        container = gammalib.GObservations()
        container.append(obs)

        # Filter event list according to ROI and user energy range
        select = ctools.ctselect(container)
        select['ra'] = ra
        select['dec'] = dec
        select['rad'] = rad
        select['emin'] = self['emin'].real()
        select['emax'] = self['emax'].real()
        select['tmin'] = 'INDEF'
        select['tmax'] = 'INDEF'
        select.run()

        # Get deep copy of filtered observation
        new_obs = select.obs()[0].copy()

        # Return
        return new_obs

    def _mask_onoff(self,obs,ra,dec,rad):
        """
        Create On/Off observation
        with On region matching the required mask

        Parameters
        ----------
        obs     : `gammalib.GCTAObservation` Input observation of type cube
        ra      : `float` R.A. (deg)
        dec     : `float` Dec (deg)
        rad     : `float` radius (deg)

        Returns
        -------
        new_obs : `gammalib.GCTAOnoffObservation`
            Output On/Off observation
        """

        # Put observation in container
        container = gammalib.GObservations()
        container.append(obs)

        onoff_obs = obsutils.get_onoff_obs(self, container, ra=ra, dec=dec, rad=rad)

        return onoff_obs[0]

    def _mask_observations(self,ra,dec,rad):
        """
        Create observations restricted to circular ROI

        Parameters
        ----------
        ra      : `float` R.A. (deg)
        dec     : `float` Dec (deg)
        rad     : `float` radius (deg)

        Returns
        -------
        new_obs : `~gammalib.GObservations`
            Observations in circular ROI
        """

        # Create output observation container
        new_obs = gammalib.GObservations()

        # Loop over input observations
        for obs in self.obs():
            if obs.classname() == 'GCTAObservation':
                if obs.eventtype() == 'CountsCube':
                    masked_obs = self._mask_cube(obs,ra,dec,rad)
                else:
                    if self['method'].string() == '3D':
                        masked_obs = self._mask_evlist(obs,ra,dec,rad)
                    elif self['method'].string() == 'ONOFF':
                        pass
                        # TODO: implement Onoff analysis
                        # masked_obs = self._mask_onoff(obs, ra, dec, rad)
                new_obs.append(masked_obs)
            # Skip On-Off and non-CTA observations
            else:
                pass

        # Append models to output observations
        new_obs.models(self.obs().models())

        # Return
        return new_obs

    def _pixel_analysis(self,ra,dec):
        """
        Performs analysis over the region of interest
        corresponding to a single pixel of output map

        Parameters
        ----------
        inx     : `int` pixel index
        ra      : `float` R.A. (deg)
        dec     : `float` Dec (deg)

        Returns
        -------
        result  : `dict` Results
        """

        # Write header for spatial bin
        msg = 'Spatial bin centred on (R.A.,Dec) = (%f,%f) deg' % (ra,dec)
        self._log_header2(gammalib.EXPLICIT, msg)

        # Initialize results
        result = {}
        for name in self._srcnames:
            result[name] = {'flux'    : 0.0,
                            'flux_err': 0.0,
                            'TS'      : 0.0,
                            'ulimit'  : 0.0}

        # Mask observations
        self._log_header3(gammalib.EXPLICIT, 'Masking observations')
        masked_obs = self._mask_observations(ra,dec,self._roisz)

        # Set up likelihood analysis
        self._log_header3(gammalib.EXPLICIT, 'Performing fit in energy bin')
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

        # Prepare objects for flux extraction
        # ROI
        centre = gammalib.GSkyDir()
        centre.radec_deg(ra,dec)
        roi = gammalib.GSkyRegionCircle(centre,self._roisz)
        # Energy boundaries
        emin = gammalib.GEnergy(self['emin'].real(),'TeV')
        emax = gammalib.GEnergy(self['emax'].real(), 'TeV')

        # Continue only if log-likelihood is non-zero
        if like.obs().logL() != 0.0:

            # Loop over target sources
            for name in self._srcnames:

                # Get source
                source = like.obs().models()[name]

                # Get flux
                # Integrate spectral model between emin and emax
                flux = source.spectral().flux(emin,emax)
                # Multiply by flux from spatial component in ROI
                flux *= source.spatial().flux(roi)
                # Divide by solid angle of ROI
                flux /= 1 - math.cos(math.radians(self._roisz))
                result[name]['flux'] = flux

                # Get flux error
                # Normalization parameter
                normpar = source.spectral()[0]
                # Only normalization parameter free
                # Relative error on flux is same as on normalization
                # Avoid zero division error
                if normpar.value() > 0.:
                    flux_error = flux * normpar.error()/normpar.value()
                else:
                    flux_error = 0.
                result[name]['flux_error'] = flux_error

                # If requested get TS
                if self['calc_ts'].boolean():
                    result[name]['TS'] = source.ts()

                # Upper limit
                # TODO

        else:
            value = 'Likelihood is zero. Bin is skipped.'
            self._log_value(gammalib.TERSE, '(R.A.,Dec) = (%f,%f) deg' % (ra,dec), value)

        # Return result
        return result

    def _fill_fits(self,results):
        """
        Fill FITS object to store the results

        Parameters
        ----------
        result  : `dict` Results

        Returns
        -------
        fits    : `~gammalib.GFits` Fits object
        """

        # Create new Fits object
        fits = gammalib.GFits()

        # If more than one source of interest
        # append empty primary HDU
        if len(self._srcnames) > 1:
            fits.append(gammalib.GFitsImageDouble())

        # Prepare some header cards
        # Flux units
        fcard = gammalib.GFitsHeaderCard('BUNIT', 'ph/cm2/s/sr' , 'Photon flux')
        txt = 'Correlation radius %f deg' %(self._roisz)
        roiszcard = gammalib.GFitsHeaderCard('COMMENT', '' , txt)

        # Loop over target sources
        for s, name in enumerate(self._srcnames):

            # Create minimal set of skymaps to store fluxes and flux errors
            fmap = self._create_map()
            errmap = self._create_map()

            #Create map for TS
            if self['calc_ts'].boolean():
                tsmap = self._create_map()

            # Loop over pixels and fill maps
            for inx in range(fmap.npix()):

                # Fill flux and flux error
                fmap[inx,0] = results[inx][name]['flux']
                errmap[inx, 0] = results[inx][name]['flux_error']

                # If requested fill TS
                if self['calc_ts'].boolean():
                    tsmap[inx, 0] = results[inx][name]['TS']

            # Write maps to Fits
            fhdu = fmap.write(fits, name + ' FLUX')
            errhdu = errmap.write(fits, name + ' FLUX ERROR')

            # Add units and comments to headers
            fhdu.header().append(fcard)
            errhdu.header().append(fcard)
            fhdu.header().append(roiszcard)
            errhdu.header().append(roiszcard)

            # If requested write TS map to fits
            if self['calc_ts'].boolean():
                tshdu = tsmap.write(fits, name + ' TS')
                tshdu.header().append(roiszcard)

        # Return
        return fits


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

        # Create template for output map, to use to extract pixel values
        outmap = self._create_map()

        # Adjust model parameters dependent on input user parameters
        self._adjust_model_pars()

        # Loop over pixels and extract results
        results = {}
        for inx in range(outmap.npix()):

            # Determine pixel centre
            pixdir = outmap.inx2dir(inx)
            ra = pixdir.ra_deg()
            dec = pixdir.dec_deg()

            # Analyse ROI
            result = self._pixel_analysis(ra,dec)
            results[inx] = result

        # Fill output Fits
        self._fits = self._fill_fits(results)

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
        # If a regions object is provided then set the regions ...
        if object is not None:
            self._excl_reg_map = gammalib.GSkyRegionMap(object)

        # Return
        return self._excl_reg_map


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csscs(sys.argv)

    # Execute application
    app.execute()
