#!/usr/bin/env python
# ==========================================================================
# Selects model from a model definition XML file
#
# Copyright (C) 2016 Juergen Knoedlseder
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


# =================== #
# csmodelselect class #
# =================== #
class csmodelselect(ctools.cscript):
    """
    Selects model from a model definition XML file
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'csmodelselect'
        self._version = '1.2.0'

        # Initialise observation container from constructor arguments
        self._obs, argv = self._set_input_obs(argv)

        # Initialise class members
        self._models = gammalib.GModels()

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Query input parameters
        self['inobs'].filename()
        self['inmodel'].filename()

        # Query hidden parameters
        self['roilimit'].real()
        self['roimargin'].real()
        self['ethres'].real()
        self['fluxlimit'].real()
        self['tslimit'].real()
        self['fit_pos'].boolean()
        self['fit_shape'].boolean()

        # Query ahead output model filename
        if self._read_ahead():
            self['outmodel'].filename()

        # If there are no observations in container then get them from the
        # parameter file
        if self._obs.size() == 0:
            self._obs = self._get_observations(False)

        # Get models
        self._models = gammalib.GModels(self['inmodel'].filename())

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _select_model(self, model, obs):
        """
        Select model

        Parameters
        ----------
        model : `~gammalib.GModel`
            Model
        obs : `~gammalib.GObservations`
            Observation container

        Returns
        -------
        select : bool
            Model selection flag
        """
        # Get selection parameters
        roilimit  = self['roilimit'].real()
        roimargin = self['roimargin'].real()
        ethres    = self['ethres'].real()
        fluxlimit = self['fluxlimit'].real()
        tslimit   = self['tslimit'].real()

        # Initialise selection flag to True
        select = True
        msg    = 'Select by default'

        # If the model has a spatial component then check of overlap
        if hasattr(model, 'spatial'):

            # Determine the model bounding box
            model_bounds = model.spatial().region()

            # Initialise with an unselect model
            select = False
            msg    = 'Exclude since outside any RoI'

            # Loop over all observations and check whether the model falls
            # into one of them
            for o in obs:

                # Get RoI centre and radius. The RoI radius is limited by
                # roilimit. A margin given by roimargin is added to the RoI
                # radius.
                #obs_centre = o.roi().centre().dir() # Segmentation fault
                #obs_radius = o.roi().radius()       # Segmentation fault
                obs_roi    = o.roi()
                obs_centre = obs_roi.centre().dir()
                obs_radius = obs_roi.radius()
                if obs_radius > roilimit:
                    obs_radius = roilimit
                obs_radius += roimargin

                # Set circular sky region
                obs_bounds = gammalib.GSkyRegionCircle(obs_centre, obs_radius)

                # If model overlaps with RoI then signal overlap
                if obs_bounds.overlaps(model_bounds):
                    select = True
                    msg    = 'Select due to overlap with at least one RoI'
                    break

        # If model is selected and if model has a TS value then apply TS
        # selection
        if select and model.has_ts():
            if model.ts() < tslimit:
                select = False
                msg    = 'Exclude since below TS limit (TS=%.3f)' % model.ts()

        # If model is selected and if model is a sky model then apply flux
        # limit selection
        if select and model.classname() == 'GModelSky':
            emin = gammalib.GEnergy(ethres,   'TeV')
            emax = gammalib.GEnergy(1000.0, 'TeV')
            flux = model.spectral().flux(emin, emax)
            if flux < fluxlimit:
                select = False
                msg    = 'Exclude since below flux limit (Flux=%.3e ph/cm2/s)' % flux

        # Log model selection
        self._log_value(gammalib.NORMAL, model.name(), msg)

        # Return selection flag
        return select

    def _set_model_parameters(self, model):
        """
        Set model parameters

        Parameters
        ----------
        model : `~gammalib.GModel`
            Model

        Returns
        -------
        model : `~gammalib.GModel`
            Model with parameter set
        """
        # Get selection parameters
        fit_pos   = self['fit_pos'].boolean()
        fit_shape = self['fit_shape'].boolean()

        # If the model has a spatial component then set the spatial parameters
        if hasattr(model, 'spatial'):

            # Loop over all parameters
            for par in model:

                # Handle position parameters
                if par.name() == 'RA' or par.name() == 'DEC':
                    if fit_pos:
                        par.free()
                    else:
                        par.fix()

                # Handle shape parameters
                elif par.name() == 'Radius' or \
                     par.name() == 'Sigma' or \
                     par.name() == 'Width' or \
                     par.name() == 'PA' or \
                     par.name() == 'MajorRadius' or \
                     par.name() == 'MinorRadius':
                    if fit_shape:
                        par.free()
                    else:
                        par.fix()

        # Return model
        return model


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

        # Initialise empty model container for selected models
        models = gammalib.GModels()

        # Write input observation container into logger
        self._log_observations(gammalib.NORMAL, self._obs, 'Input observation')

        # Write input models into logger
        self._log_models(gammalib.VERBOSE, self._models, 'Input model')

        # Write header
        self._log_header1(gammalib.NORMAL, 'Model selection')

        # Loop over model components
        for model in self._models:

            # If model should be selected then set model parameters and append
            # model to the output container
            if self._select_model(model, self._obs):

                # Set model parameters
                model = self._set_model_parameters(model)

                # Append model
                models.append(model)

        # Copy selected model into selected models
        self._models = models

        # Write selected models into logger
        self._log_models(gammalib.VERBOSE, self._models, 'Selected model')

        # Return
        return

    def save(self):
        """ 
        Save model definition XML file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save models')

        # Get output filename in case it was not read ahead
        outmodel = self['outmodel'].filename()

        # If file exists and clobber flag is false then raise an exception
        if outmodel.exists() and not self._clobber:
            msg = ('Cannot save "'+outmodel.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:
            # Log filename
            self._log_value(gammalib.NORMAL, 'Model definition XML file',
                                             outmodel.url())

            # Save models
            self._models.save(outmodel)

        # Return
        return

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save model file if required
        self.save()

        # Return
        return    


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csmodelselect(sys.argv)

    # Execute application
    app.execute()
