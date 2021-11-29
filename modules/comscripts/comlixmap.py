#! /usr/bin/env python
# ==========================================================================
# Create SRCLIX TS map
#
# Copyright (C) 2021 Juergen Knoedlseder
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
        self._maps      = []
        self._map_names = []
        self._srcnames  = []

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

        # Query parameters
        self['like_accuracy'].real()
        self['max_iter'].integer()
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

        # Compute number of free model parameters for all test sources,
        # exclusing 'RA' and 'DEC'
        for srcname in self._srcnames:
            source = self.obs().models()[srcname]
            for par in source:
                if par.name() != 'RA' and par.name() != 'DEC' and par.is_free():
                    self._map_names.append(par.name())
                    nfree += 1

        # Initialise sky maps
        for i in range(nfree+len(self._srcnames)):
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
        # Loop over test source names
        for srcname in self._srcnames:

            # Set and fix position of test source
            source = self.obs().models()[srcname]
            source['RA'].value(dir.ra_deg())
            source['DEC'].value(dir.dec_deg())
            source['RA'].fix()
            source['DEC'].fix()

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

            # Fit model
            self.obs().optimize(self.opt())

            # Compute logL difference after first iteration
            if iter > 0:
                delta = logL - self.opt().value()

            # Store maximum likelihood value
            logL = self.opt().value()

            # Log maximum likelihood
            if iter == 0:
                result = '%.5f' % (logL)
            else:
                result = '%.5f (%.5f)' % (logL, delta)
            self._log_value(gammalib.NORMAL, 'logL after iteration %d' % (iter+1),
                            result)

            # Check for convergence
            if iter > 0:
                if delta < eps:
                    break

        # Compute TS in final model fit
        ts = self._final_model_fit()

        # Log TS
        for i, srcname in enumerate(self._srcnames):
            key = 'TS %s' % srcname
            self._log_value(gammalib.NORMAL, key, ts[i])

        # Return TS
        return ts

    def _final_model_fit(self):
        """
        Perform final model fit using ctlike

        Returns
        -------
        ts : float
            TS of test source
        """
        # Create instance of model fitting tool
        like = ctools.ctlike(self.obs())
        like['fix_spat_for_ts'] = self['fix_spat'].boolean()

        # Run ctlike
        like.run()

        # Recover TS for all test sources
        ts = []
        for srcname in self._srcnames:
            ts.append(like.obs().models()[srcname].ts())

        # Return TS
        return ts


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
        self._log_header1(gammalib.NORMAL, 'Initialise TS map')

        # Create sky maps
        self._create_maps()

        # Log sky map
        self._log_string(gammalib.NORMAL, str(self._maps[0]))

        # Log header
        self._log_header1(gammalib.NORMAL, 'Generate TS map')

        # Loop over grid positions
        for ipix in range(self._maps[0].npix()):

            # Get sky direction of sky map pixel
            dir = self._maps[0].inx2dir(ipix)

            # Log header
            header = 'Computing TS for pixel %d/%d at %s' % \
                     (ipix, self._maps[0].npix(), str(dir))
            self._log_header3(gammalib.NORMAL, header)

            # Compute TS
            ts = self._compute_ts(dir)

            # Store TS values in sky map
            for i, value in enumerate(ts):
                self._maps[i][ipix] = value

            # Store fitted model parameters in sky maps
            ipar   = len(ts)
            for srcname in self._srcnames:
                source = self.obs().models()[srcname]
                for par in source:
                    if par.name() != 'RA' and par.name() != 'DEC' and par.is_free():
                        self._maps[ipar][ipix] = par.value()
                        ipar += 1

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
            for i, srcname in enumerate(self._srcnames):
                fits[i].extname(srcname)
            for i, name in enumerate(self._map_names):
                fits[i+len(self._srcnames)].extname(name)

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
