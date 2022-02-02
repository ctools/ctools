#! /usr/bin/env python
# ==========================================================================
# Perform SRCLIX model fitting of COMPTEL observations
#
# Copyright (C) 2021-2022 Juergen Knoedlseder
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
import os
import sys
import gammalib
import ctools


# =============== #
# comlixfit class #
# =============== #
class comlixfit(ctools.cslikelihood):
    """
    Perform SRCLIX model fitting of COMPTEL observations
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_cslikelihood(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._fixed_pars = []
        self._iter       = 0
        self._logL       = 0.0
        self._nobs       = 0.0
        self._npred      = 0.0

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

        # Query parameters
        self['max_iter'].integer()
        self['like_accuracy'].real()
        self['fix_spat_for_ts'].boolean()
        self['fix_ext_for_initial'].boolean()

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
            self['suffix'].string()
            self['outfolder'].string()
            self['outobs'].filename()
            self['outmodel'].filename()

        # Set Phibar layers for fitting
        for obs in self.obs():
            if obs.classname() == 'GCOMObservation':
                if phi_first >= 0:
                    obs.phi_first(phi_first)
                if phi_last >= 0:
                    obs.phi_last(phi_last)

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Set optimizer logger
        if self._logExplicit():
            self.opt().logger(self._log)

        # Return
        return

    def _update_obs(self):
        """
        Update background model in observation container

        The method updated the background model in the observation container
        by taking into account the current source models in the BGDLIXA
        model generation algorithm.
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

            # Skip non-COMPTEL observations
            if obs.classname() != 'GCOMObservation':
                continue

            # Compute DRM
            drm = obs.drm(models)

            # Compute background model
            obs.compute_drb(bkgmethod, drm, nrunav, navgr, nincl, nexcl)

            # Signal that DRB file was not yet saved
            obs.drbname('')

        # Return
        return

    def _final_model_fit(self):
        """
        Perform final model fit using ctlike
        """
        # Create instance of model fitting tool
        like = ctools.ctlike(self.obs())
        like['fix_spat_for_ts'] = self['fix_spat_for_ts'].boolean()

        # Run ctlike
        like.run()

        # Recover results
        self.opt(like.opt())
        self.obs(like.obs())
        self._iter  += like.iter()
        self._logL   = like.logL()
        self._nobs   = like.nobs()
        self._npred  = like.npred()

        # Return
        return

    def _fix_ext_pars(self):
        """
        Fix spatial extension parameters
        """
        # Initialise list of fixed parameters
        self._fixed_pars = []

        # Set list of spatial extension parameters
        pars = ['Radius', 'Width', 'Sigma', 'PA', 'MajorRadius', 'MinorRadius']

        # Loop over all sky models
        for model in self.obs().models():
            if model.classname() == 'GModelSky':

                # Loop over all parameters and fix all encountered spatial
                # extension parameters
                for par in pars:
                    if model.spatial().has_par(par):
                        if model[par].is_free():

                            # Fix parameter
                            model[par].fix()

                            # Bookkeeping fixed parameters
                            fixed = {'name': model.name(), 'par':  par}
                            self._fixed_pars.append(fixed)

                            # Log action
                            value = '%s[%s]' % (model.name(), par)
                            self._log_value(gammalib.EXPLICIT, 'Fixed parameter', value)

        # Return
        return

    def _free_fixed_ext_pars(self):
        """
        Free fixed spatial extension parameters
        """
        # Continue only if there are fixed spatial extension parameters
        if len(self._fixed_pars) > 0:

            # Loop over fixed parameters
            for fixed_par in self._fixed_pars:

                # Free parameter
                name = fixed_par['name']
                par  = fixed_par['par']
                self.obs().models()[name][par].free()

                # Log action
                value = '%s[%s]' % (name, par)
                self._log_value(gammalib.EXPLICIT, 'Free fixed parameter', value)

            # Clear list of fixed parameters
            self._fixed_pars = []

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

        # Log header
        self._log_header1(gammalib.NORMAL, 'Input observations')

        # Log input observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

        # Get parameters and initialise some variables
        niter      = self['max_iter'].integer()
        eps        = self['like_accuracy'].real()
        delta      = 0.0
        self._iter = 0

        # Log header
        self._log_header1(gammalib.NORMAL,
                          'Iterative maximum likelihood model fitting')

        # Loop over iterations
        for iter in range(niter):

            # Log iteration header
            self._log_header2(gammalib.EXPLICIT, 'Iteration %d' % (iter+1))

            # Update observations
            self._update_obs()

            # Optionally fix or free fixed spatial extension parameters
            if self['fix_ext_for_initial'].boolean():
                if iter == 0:
                    self._fix_ext_pars()
                else:
                    self._free_fixed_ext_pars()

            # Fit model
            self.obs().optimize(self.opt())

            # Compute logL difference after first iteration
            if iter > 0:
                delta = logL - self.opt().value()

            # Store maximum likelihood value
            logL = self.opt().value()

            # Increment iterations counter
            self._iter += self.opt().iter()

            # Log maximum likelihood
            if iter == 0:
                result = '%.5f' % (logL)
            else:
                result = '%.5f (%.5f)' % (logL, delta)
            self._log_value(gammalib.NORMAL, 'logL after iteration %d' % (iter+1),
                            result)

            # Optionally log results
            self._log_string(gammalib.VERBOSE, str(self.obs().models()))

            # Check for convergence
            if iter > 0:
                if delta < eps:
                    break

        # Log header
        self._log_header1(gammalib.EXPLICIT,
                          'Final maximum likelihood model fitting')

        # Do final model fit
        self._final_model_fit()

        # Compute logL difference and store maximum likelihood value
        delta  = logL - self.opt().value()
        logL   = self.opt().value()

        # Log final maximum likelihood
        result = '%.5f (%.5f)' % (logL, delta)
        self._log_value(gammalib.NORMAL, 'logL after final iteration',
                        result)

        # Log results
        self._log_header1(gammalib.NORMAL,
                          'Maximum likelihood optimisation results')
        self._log_string(gammalib.NORMAL, str(self.opt()))
        self._log_value(gammalib.NORMAL, 'Total number of iterations',
                        self._iter)
        self._log_value(gammalib.NORMAL, 'Maximum log likelihood',
                        self._logL)
        self._log_value(gammalib.NORMAL, 'Observed events  (Nobs)',
                        self._nobs)
        self._log_value(gammalib.NORMAL, 'Predicted events (Npred)',
                        self._npred)
        self._log_value(gammalib.NORMAL, 'Nobs - Npred',
                        self._nobs-self._npred)
        self._log_string(gammalib.NORMAL, str(self.obs().models()))

        # Return
        return

    def save(self):
        """ 
        Save observation definition file
        """
        # Log header
        self._log_header1(gammalib.TERSE, 'Save observations')

        # Get output filenames
        outobs   = self['outobs'].filename()
        outmodel = self['outmodel'].filename()

        # If file exists and clobber flag is false then raise an exception
        if outobs.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outobs.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)
        elif outmodel.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outmodel.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:

            # Get DRB file suffix and set outfolder
            suffix    = self['suffix'].string()
            outfolder = self['outfolder'].string()

            # Create outfolder directory
            try:
                os.makedirs(gammalib.expand_env(outfolder))
            except OSError:
                pass

            # Loop over all observations
            for obs in self.obs():

                # Skip non-COMPTEL observations
                if obs.classname() != 'GCOMObservation':
                    continue

                # Store background filename
                drename = '%s/%s' % (outfolder, os.path.basename(obs.drename().url()))
                if suffix == '':
                    drbname = drename.replace('dre', 'drb')
                else:
                    drbname = drename.replace('dre', 'drb-%s' % suffix)
                obs.drbname(drbname)

                # Save DRB file
                obs.drb().save(drbname, self['clobber'].boolean())

                # Log saving
                self._log_value(gammalib.NORMAL, 'DRB file', drbname)

            # Log observation definition filename
            self._log_value(gammalib.NORMAL, 'Obs. definition XML file',
                                             outobs.url())

            # Save observations
            self.obs().save(outobs)

            # Log model definition filename
            self._log_value(gammalib.NORMAL, 'Model definition XML file',
                                             outmodel.url())

            # Save models
            self.obs().models().save(outmodel)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comlixfit(sys.argv)

    # Execute application
    app.execute()
