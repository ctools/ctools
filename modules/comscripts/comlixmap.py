#! /usr/bin/env python
# ==========================================================================
# Create SRCLIX TS map
#
# Copyright (C) 2021-2023 Juergen Knoedlseder
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


# =============== #
# comlixmap class #
# =============== #
class comlixmap(ctools.cslikelihood):
    """
    Create SRCLIX TS map
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_cslikelihood(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._inmaps      = []
        self._inmap_names = []
        self._maps        = []
        self._map_names   = []
        self._srcnames    = []

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Set observation if not done before
        if self.obs().is_empty():
            self.obs().load(self['inobs'].filename())

        # Set models if we have none
        if self.obs().models().is_empty():
            self.obs().models(self['inmodel'].filename())

        # Get test source names and check if test sources exists in model
        self._srcnames = self['srcname'].string().split(';')
        for srcname in self._srcnames:
            if not self.obs().models().contains(srcname):
                msg = 'Source "%s" not found in models.' % srcname
                raise RuntimeError(msg)

        # Query optional input TS map filename
        if self['inmap'].is_valid():
            self['inmap'].filename()

        # Query parameters
        self['like_accuracy'].real()
        self['max_iter'].integer()
        self['accept_dec'].real()
        self['fix_spat'].boolean()
        self['coordsys'].string()
        self['proj'].string()
        self['xref'].real()
        self['yref'].real()
        self['binsz'].real()
        self['nxpix'].integer()
        self['nypix'].integer()

        # Get parameters
        bkgmethod = self['bkgmethod'].string()
        nrunav    = self['nrunav'].integer()
        navgr     = self['navgr'].integer()
        nincl     = self['nincl'].integer()
        nexcl     = self['nexcl'].integer()
        phi_first = self['phi_first'].integer()
        phi_last  = self['phi_last'].integer()

        # Check for incorrect parameters
        if nexcl < 0 or nexcl >= nincl:
            msg = 'Incorrect value %d for nexcl (bins to exclude).' % nexcl
            raise RuntimeError(msg)
        if nexcl != 0 and 2*int(nexcl/2) == nexcl :
            msg = 'nexcl=%d (bins to exclude) should be zero or odd number.' % nexcl
            raise RuntimeError(msg)
        if nincl < 3 or 2*int(nincl/2) == nincl:
            msg = 'nincl=%d (bins to include) should be odd and >= 3.' % nincl
            raise RuntimeError(msg)
        if navgr < 1 or 2*int(navgr/2) == navgr :
            msg = 'navgr=%d should be odd and >= 1.' % navgr
            raise RuntimeError(msg)

        # Query ahead output model filename
        if self._read_ahead():
            self['outmap'].filename()

        # Set Phibar layers for fitting
        for obs in self.obs():
            if obs.classname() == 'GCOMObservation':
                if phi_first >= 0:
                    obs.phi_first(phi_first)
                if phi_last >= 0:
                    obs.phi_last(phi_last)

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _update_obs(self):
        """
        Update background model in observation container

        The method updates the background model in the observation container
        by taking into account the current source models in the  model
        generation algorithm.
        """
        # Get task parameters
        bkgmethod = self['bkgmethod'].string()
        nrunav    = self['nrunav'].integer()
        navgr     = self['navgr'].integer()
        nincl     = self['nincl'].integer()
        nexcl     = self['nexcl'].integer()

        # Extract source models from observation model container
        models = gammalib.GModels()
        for model in self.obs().models():
            if model.classname() == 'GModelSky' or model.classname() == 'GCOMModelDRM':
                models.append(model)

        # Loop over all observations
        for obs in self.obs():

            # Compute DRM
            drm = obs.drm(models)

            # Compute background model
            obs.compute_drb(bkgmethod, drm, nrunav, navgr, nincl, nexcl)

        # Return
        return

    def _load_inmaps(self, filename):
        """
        Load input TS maps

        Loads input TS maps from a FITS files and checks that number
        of extensions and extension names are consistent with the
        expectations.

        Parameters
        ----------
        filename : `~gammalib/GFilename`
            Input filename
        """
        # Log header
        self._log_header1(gammalib.NORMAL, 'Load input TS map(s)')
        self._log_value(gammalib.NORMAL, 'Input filename', filename.url())

        # Initialise list of input maps and map names
        self._inmaps      = []
        self._inmap_names = []

        # Open FITS file
        fits = gammalib.GFits(filename)

        # Raise an exception if the number of extensions does not
        # match the expectations
        if fits.size() != len(self._map_names):
            msg = 'There are %d extensions in the input TS map but %d '\
                  'extensions are expected. Please provide a compatible '\
                  'input TS map.' % (fits.size(), len(self._map_names))
            raise RuntimeError(msg)

        # Loop over extensions and extract all images
        for i, hdu in enumerate(fits):

            # Get extension name
            name = hdu.extname()

            # Raise exception if the extension name is conform to
            # expectation
            if name != self._map_names[i]:
                msg = 'Encountered extension name "%s" of map %d '\
                      'does not correspond to expected extension '\
                      'name "%s". Please provide a compatible '\
                      'input TS map.' % (name, i, self._map_names[i])
                raise RuntimeError(msg)

            # Append extension name and map
            map = gammalib.GSkyMap(hdu)
            self._inmap_names.append(name)
            self._inmaps.append(map)
            self._log_value(gammalib.NORMAL, 'Loaded', name)

        # Return
        return

    def _prefit_models(self):
        """
        Prefits models

        The method prefits models to the data.
        """
        # Get parameters and initialise some variables
        niter = self['max_iter'].integer()
        eps   = self['like_accuracy'].real()
        delta = 0.0

        # Loop over iterations
        for iter in range(niter):

            # Update observations
            self._update_obs()

            # If first iteration then initialise best observation
            if iter == 0:
                best_obs = self.obs().copy()

            # Fit model
            self.obs().optimize(self.opt())

            # Compute log-likelihood difference or initialise
            # log-likelihood value
            if iter > 0:
                delta = logL - self.opt().value()
            else:
                logL = self.opt().value()

            # If log-likelihood improved then update best observations and
            # last log-likelihood value
            if delta >= 0.0:
                best_obs = self.obs().copy()
                logL     = self.opt().value()

            # Log maximum likelihood
            if iter == 0:
                result = '%.5f' % (self.opt().value())
            else:
                result = '%.5f (%.5f)' % (self.opt().value(), delta)
            self._log_value(gammalib.NORMAL, 'logL after iteration %d' % (iter+1),
                            result)

            # Check for convergence
            if iter > 0:
                if delta < eps:
                    break

        # Use best observations for maximum likelihood fitting. This avoids
        # that observations that have a worse log-likelihood than the best
        # value are used for final model fitting.
        self.obs(best_obs)

        # Compute log-likelihood and TS for best observations
        value, opt = self._final_model_fit()

        # Compute log-likelihood difference
        delta = logL - value

        # Log final maximum likelihood
        result = '%.5f (%.5f)' % (value, delta)
        self._log_value(gammalib.NORMAL, 'logL after final iteration',
                        result)

        # Log optimiser
        self._log_string(gammalib.NORMAL, str(opt))

        # Return
        return


    def _create_maps(self):
        """
        Create sky maps based on user parameters

        The method initialises all sky maps for the script and stores the
        FITS extension names for the test source parameters.
        """
        # Initialise list of maps and map names
        self._maps      = []
        self._map_names = []

        # Get user parameters
        coordsys = self['coordsys'].string()
        proj     = self['proj'].string()
        xref     = self['xref'].real()
        yref     = self['yref'].real()
        binsz    = self['binsz'].real()
        nxpix    = self['nxpix'].integer()
        nypix    = self['nypix'].integer()

        # Initialise number of free parameters
        nfree = 0

        # Append TS map names
        for srcname in self._srcnames:
            self._map_names.append(srcname)

        # Compute number of free model parameters for all test sources,
        # excluding 'RA', 'DEC', 'GLON' and 'GLAT' and add map names
        # composed of 'srcname_parname_val' and 'srcname_parname_unc'
        # for values and uncertainties
        for srcname in self._srcnames:
            source = self.obs().models()[srcname]
            for par in source:
                if par.name() != 'RA'   and par.name() != 'DEC'  and \
                   par.name() != 'GLON' and par.name() != 'GLAT' and \
                   par.is_free():
                    self._map_names.append('%s_%s_val' % (srcname, par.name()))
                    self._map_names.append('%s_%s_unc' % (srcname, par.name()))
                    nfree += 1

        # Initialise sky maps
        for i in range(2*nfree+len(self._srcnames)):
            self._maps.append(gammalib.GSkyMap(proj, coordsys, xref, yref,
                                               -binsz, binsz, nxpix, nypix))

        # Return
        return

    def _compute_ts(self, dir):
        """
        Compute TS

        The method computes the TS value for a given test source sky direction.

        Parameters
        ----------
        dir : `~gammalib.GSkyDir`
            Position of test source

        Returns
        -------
        ts : float
            TS of test source
        """
        # Get list of flux normalisation parameters
        pars = ['Prefactor', 'Normalization', 'EnergyFlux', 'PhotonFlux']
        
        # Loop over test source names
        for srcname in self._srcnames:

            # Set and fix position of test source
            source = self.obs().models()[srcname]
            if source.has_par('RA') and source.has_par('DEC'):
                source['RA'].value(dir.ra_deg())
                source['DEC'].value(dir.dec_deg())
                source['RA'].fix()
                source['DEC'].fix()
            elif source.has_par('GLON') and source.has_par('GLAT'):
                source['GLON'].value(dir.l_deg())
                source['GLAT'].value(dir.b_deg())
                source['GLON'].fix()
                source['GLAT'].fix()
            else:
                msg = ('Source "'+srcname+'" has neither RA/DEC nor GLON/GLAT '
                       'parameters. Please specify a source than can be '
                       'positioned.')
                raise RuntimeError(msg)

            # Remove response cache for test source
            self.obs().remove_response_cache(srcname)

        # Get parameters and initialise some variables
        niter = self['max_iter'].integer()
        eps   = self['like_accuracy'].real()
        delta = 0.0

        # Loop over iterations
        for iter in range(niter):

            # Update observations
            self._update_obs()

            # If first iteration then initialise best observation
            if iter == 0:
                best_obs = self.obs().copy()

            # Fit model
            self.obs().optimize(self.opt())

            # Compute log-likelihood difference or initialise
            # log-likelihood value
            if iter > 0:
                delta = logL - self.opt().value()
            else:
                logL = self.opt().value()

            # If log-likelihood improved then update best observations and
            # last log-likelihood value
            if delta >= 0.0:
                best_obs = self.obs().copy()
                logL     = self.opt().value()

            # Log maximum likelihood
            if iter == 0:
                result = '%.5f' % (self.opt().value())
            else:
                result = '%.5f (%.5f)' % (self.opt().value(), delta)
            self._log_value(gammalib.NORMAL, 'logL after iteration %d' % (iter+1),
                            result)

            # Check for convergence
            if iter > 0:
                if delta < eps:
                    break

        # Use best observations for maximum likelihood fitting. This avoids
        # that observations that have a worse log-likelihood than the best
        # value are used for final model fitting.
        self.obs(best_obs)

        # Compute log-likelihood and TS for best observations
        value, _ = self._final_model_fit()

        # Compute log-likelihood difference
        delta = logL - value

        # Log final maximum likelihood
        result = '%.5f (%.5f)' % (value, delta)
        self._log_value(gammalib.NORMAL, 'logL after final iteration',
                        result)

        # Collect and log TS and source normalisation
        ts = []
        for srcname in self._srcnames:
            source   = self.obs().models()[srcname]
            ts_value = float(source.ts())
            key      = 'TS %s' % srcname
            value    = '%.3f'  % ts_value
            ts.append(ts_value)
            for par in pars:
                if source.has_par(par):
                    value += ' (%e +/- %e %s)' % (source[par].value(),
                                                  source[par].error(),
                                                  source[par].unit())
                    break
            self._log_value(gammalib.NORMAL, key, value)

        # Return TS
        return ts

    def _extract_from_inmaps(self, dir, ipix):
        """
        Extract information from input maps

        Parameters
        ----------
        dir : `~gammalib.GSkyDir`
            Position of test source
        ipix : int
            TS map pixel

        Returns
        -------
        extracted : boolean
            True if extraction was successful
        """
        # Initialise extraction flag
        extracted = False

        # If there are input TS maps then check if the direction is
        # contained in the maps
        if len(self._inmaps) == len(self._maps):

            # Get first input TS map
            map = self._inmaps[0]

            # If position of test source is contained in input TS map
            # then check whether the corresponding pixel distance is
            # sufficiently small
            if map.contains(dir):

                # If distance between the position of the test source
                # to input TS map pixel is smaller than 0.1 degrees
                # then recover information
                inx     = map.dir2inx(dir)
                map_dir = map.inx2dir(inx)
                dist    = map_dir.dist_deg(dir)
                if dist < 0.1:

                    # Recover information
                    for i in range(len(self._inmaps)):
                        self._maps[i][ipix] = self._inmaps[i][inx]

                    # Log that information was extraced
                    key   = 'Input map pixel'
                    value = '%d (offset: %f deg)' % (inx, dist)
                    self._log_value(gammalib.NORMAL, key, value)
                    for i in range(len(self._maps)):
                        if i < len(self._srcnames):
                            key   = 'TS %s' % (self._map_names[i])
                            value = '%.3f'  % (self._maps[i][ipix])
                        else:
                            key   = '%s' % (self._map_names[i])
                            value = '%e' % (self._maps[i][ipix])
                        self._log_value(gammalib.NORMAL, key, value)

                    # Signal that information was extracted
                    extracted = True

        # Return
        return extracted

    def _final_model_fit(self):
        """
        Perform final model fit using ctlike

        Returns
        -------
        logL : float
            Log-likelihood of final model fit
        opt : `~gammalib.GOptimizer`
            Optimizer
        """
        # Create instance of model fitting tool
        like = ctools.ctlike(self.obs())
        like['fix_spat_for_ts'] = self['fix_spat'].boolean()

        # Run ctlike
        like.run()

        # Recover results
        self.obs(like.obs())

        # Return log-likelihood value and optimiser
        return (like.opt().value(), like.opt().copy())


    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Log header
        self._log_header1(gammalib.NORMAL, 'Input observations')

        # Log input observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

        # Log header
        self._log_header1(gammalib.NORMAL, 'Initialise models')

        # Make sure that tscalc flag is only set for test sources
        for source in self.obs().models():
            if source.name() in self._srcnames:
                source.tscalc(True)
            else:
                source.tscalc(False)

        # Fix spatial parameters if requested
        if self['fix_spat'].boolean():
            for source in self.obs().models():
                if source.classname() == 'GModelSky':
                    for par in source.spatial():
                        par.fix()

        # Log models
        self._log_string(gammalib.NORMAL, str(self.obs().models()))

        # Log header
        self._log_header1(gammalib.NORMAL, 'Prefit models without test source')

        # Store initial models
        models_orig = self.obs().models().copy()

        # Remove test sources from models
        models = self.obs().models().copy()
        for srcname in self._srcnames:
            models.remove(srcname)
        self.obs().models(models)

        # Prefit models
        self._prefit_models()

        # Log prefitted models
        self._log_string(gammalib.NORMAL, str(self.obs().models()))

        # Append removed test sources to prefitted model
        models = self.obs().models().copy()
        for srcname in self._srcnames:
            models.append(models_orig[srcname])
        self.obs().models(models)

        # Log header
        self._log_header1(gammalib.NORMAL, 'Initialise TS map')

        # Create sky maps
        self._create_maps()

        # Log sky map
        self._log_string(gammalib.NORMAL, str(self._maps[0]))

        # Optionally load input TS map(s)
        if self['inmap'].is_valid():
            self._load_inmaps(self['inmap'].filename())
            self._log_string(gammalib.NORMAL, str(self._inmaps[0]))

        # Log header
        self._log_header1(gammalib.NORMAL, 'Generate TS map')

        # Get a copy of the initial models
        models = self.obs().models().copy()

        # Initialise optimiser
        self.opt().max_iter(100)
        self.opt().eps(0.005)
        self.opt().accept_dec(self['accept_dec'].real())

        # Loop over grid positions
        for ipix in range(self._maps[0].npix()):

            # Set initial models so that each pixel starts from the same
            # initial model
            self.obs().models(models)

            # Get sky direction of sky map pixel
            dir = self._maps[0].inx2dir(ipix)

            # Log header
            header = 'Computing TS for pixel %d/%d at %s' % \
                     (ipix, self._maps[0].npix(), str(dir))
            self._log_header3(gammalib.NORMAL, header)

            # Compute TS if information could not be recovered from input
            # maps
            if not self._extract_from_inmaps(dir, ipix):

                # Compute TS
                ts = self._compute_ts(dir)

                # Store TS values in sky map
                for i, value in enumerate(ts):
                    self._maps[i][ipix] = value

                # Store fitted model parameters and uncertainties in sky maps
                ipar = len(ts)
                for srcname in self._srcnames:
                    source = self.obs().models()[srcname]
                    for par in source:
                        if par.name() != 'RA'   and par.name() != 'DEC'  and \
                           par.name() != 'GLON' and par.name() != 'GLAT' and \
                            par.is_free():
                            self._maps[ipar][ipix]   = par.value()
                            self._maps[ipar+1][ipix] = par.error()
                            ipar += 2

        # Return
        return

    def save(self):
        """ 
        Save observation definition file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save TS map')

        # Get output filenames
        outmap = self['outmap'].filename()

        # If file exists and clobber flag is false then raise an exception
        if outmap.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outmap.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:

            # Log observation definition filename
            self._log_value(gammalib.NORMAL, 'TS map file', outmap.url())

            # Create fits file
            fits = gammalib.GFits()

            # Write the sky maps to the FITS file
            for map in self._maps:
                map.write(fits)

            # Set extension name for all maps
            for i, name in enumerate(self._map_names):
                fits[i].extname(name)

            # Save FITS file
            fits.saveto(outmap, self['clobber'].boolean())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comlixmap(sys.argv)

    # Execute application
    app.execute()
