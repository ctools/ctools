#! /usr/bin/env python
# ==========================================================================
# Generates background model
#
# Copyright (C) 2018 Juergen Knoedlseder
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


# ================ #
# csbkgmodel class #
# ================ #
class csbkgmodel(ctools.csobservation):
    """
    Generates background model
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Set members
        self._models     = gammalib.GModels()
        self._instrument = ''

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Setup observations (do not require response, allow event list and
        # counts cube)
        self._setup_observations(self.obs(), False, True, True)

        # Set instrument name
        self._instrument = self._get_instrument()

        # Query input parameters
        spatial = self['spatial'].string()
        if self['spatial'].string() == 'GAUSS':
            self['gradient'].boolean()
        self['spectral'].string()
        self['runwise'].boolean()
        self['rad'].real()

        # Query ahead output model filename
        if self._read_ahead():
            self['outmodel'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _get_instrument(self):
        """
        Get instrument name

        Extracts the instrument name from the observation container. If there
        are multiple instruments in the container then use the "instrument"
        parameter to set the instrument name.

        Returns
        -------
        instrument : str
            Instrument name
        """
        # Initialise instrument
        instrument = ''

        # Loop over all observations
        for obs in self.obs():
            if instrument == '':
                instrument = obs.instrument()
            elif instrument != obs.instrument():
                instrument = self['instrument'].string()
                break

        # Return instrument
        return instrument

    def _generate_initial_model(self):
        """
        Generate initial background model

        Returns
        -------
        model : `~gammalib.GModelData()`
            Background model
        """
        # Handle IRF model
        if self['spatial'].string() == 'IRF':
            epivot   = gammalib.GEnergy(1.0, 'TeV')
            spectral = gammalib.GModelSpectralPlaw(1.0, 0.0, epivot)
            model    = gammalib.GCTAModelIrfBackground(spectral)

        # Handle AEFF model
        elif self['spatial'].string() == 'AEFF':
            epivot   = gammalib.GEnergy(1.0, 'TeV')
            spectral = gammalib.GModelSpectralPlaw(1.0e-13, -2.5, epivot)
            model    = gammalib.GCTAModelAeffBackground(spectral)

        # Handle GAUSS model
        elif self['spatial'].string() == 'GAUSS':
            if self['gradient'].boolean():
                spatial = gammalib.GCTAModelSpatialMultiplicative()
                factor1 = gammalib.GCTAModelRadialGauss(3.0)
                factor2 = gammalib.GCTAModelSpatialGradient()
                spatial.append(factor1)
                spatial.append(factor2)
            else:
                spatial = gammalib.GCTAModelRadialGauss(3.0)
            epivot   = gammalib.GEnergy(1.0, 'TeV')
            spectral = gammalib.GModelSpectralPlaw(3.0e-4, -1.5, epivot)
            spectral['Prefactor'].min(1.0e-8)
            model    = gammalib.GCTAModelBackground(spatial, spectral)

        # Any other strings (should never occur)
        else:
            model = None

        # Set background name and instrument
        if model is not None:
            model.name('Background')
            model.instruments(self._instrument)

        # Return model
        return model

    def _select_events(self, obs):
        """
        Select events within a given RoI radius

        Parameters
        ----------
        obs : `~gammalib.GObservations()`
            Observation container

        Returns
        -------
        obs : `~gammalib.GObservations()`
            Observation container
        """
        # Setup task parameters
        select = ctools.ctselect(obs)
        select['ra']       = 'UNDEF'
        select['dec']      = 'UNDEF'
        select['rad']      = self['rad'].real()
        select['tmin']     = 'UNDEF'
        select['tmax']     = 'UNDEF'
        select['emin']     = 'UNDEF'
        select['emax']     = 'UNDEF'
        select['usethres'] = 'DEFAULT'

        # Select events
        select.run()

        # Extract observation
        obs = select.obs().copy()

        # Return observation
        return obs

    def _create_nodes(self, model):
        """
        Replace spectral model by node function

        Parameters
        ----------
        model : `~gammalib.GModelData()`
            Input background model

        Returns
        -------
        model : `~gammalib.GModelData()`
            Background model with spectral node function
        """
        # Extract spectral model of background component
        spectrum = model.spectral()

        # Define 8 node energies
        energies = gammalib.GEnergies()
        energies.append(gammalib.GEnergy(0.1, 'TeV'))
        energies.append(gammalib.GEnergy(0.2, 'TeV'))
        energies.append(gammalib.GEnergy(0.4, 'TeV'))
        energies.append(gammalib.GEnergy(0.8, 'TeV'))
        energies.append(gammalib.GEnergy(1.6, 'TeV'))
        energies.append(gammalib.GEnergy(3.2, 'TeV'))
        energies.append(gammalib.GEnergy(6.4, 'TeV'))
        energies.append(gammalib.GEnergy(12.8, 'TeV'))

        # Create node spectrum
        nodes = gammalib.GModelSpectralNodes(spectrum, energies)
        nodes.autoscale()

        # Set minimum node intensities (this avoids NaNs during model fitting)
        for i in range(nodes.size()):
            nodes[i].min(1.0e-6*nodes[i].value())

        # Set node spectrum
        model.spectral(nodes)

        # Return
        return model

    def _generate_bkg(self, obs):
        """
        Generate background models

        Parameters
        ----------
        obs : `~gammalib.GObservations()`
            Observations container

        Returns
        -------
        model : `~gammalib.GModelData()`
            Background model component
        """
        # Write header for event selection
        self._log_header3(gammalib.EXPLICIT, 'Select events from observation')

        # Select events
        obs = self._select_events(obs)

        # Write header for initial background model generation
        self._log_header3(gammalib.EXPLICIT, 'Generate initial background model')

        # Generate initial background model
        model = self._generate_initial_model()

        # Attach initial background model
        models = gammalib.GModels()
        models.append(model)
        obs.models(models)

        # Write header for initial model fitting
        self._log_header3(gammalib.EXPLICIT, 'Fit initial background model')

        # Perform maximum likelihood fitting with initial model
        like = ctools.ctlike(obs)
        like.run()

        # Extract fitted model
        model = like.obs().models()[0].copy()

        # If a HESS model is requested then refit a node spectrum
        if self['spectral'].string() == 'HESS':
        
            # Create nodes spectrum from fitted initial model
            model = self._create_nodes(model)

            # Attach node spectrum
            models = gammalib.GModels()
            models.append(model)
            obs.models(models)

            # Write header for node model fitting
            self._log_header3(gammalib.EXPLICIT, 'Fit nodes background model')

            # Perform maximum likelihood fitting with node model
            like = ctools.ctlike(obs)
            like.run()

            # Extract fitted model
            model = like.obs().models()[0].copy()

        # Return model
        return model

    def _generate_runwise_bkg(self):
        """
        Generate background models

        Parameters
        ----------
        responses : list of dict
            List of response dictionaries
        energy : float
            Energy to insert (TeV)
        comment : str
            Reason for energy insertion
        """
        # Loop over observations
        for run in self.obs():

            # Write header
            self._log_header2(gammalib.TERSE, self._get_obs_header(run))

            # Build observation container with single run
            obs = gammalib.GObservations()
            obs.append(run)

            # Generate background model for
            model = self._generate_bkg(obs)

            # Set model attributes
            model.name(('Background_%s' % run.id()))
            model.ids(run.id())

            # Write model
            self._log_string(gammalib.NORMAL, str(model))

            # Append background model
            self._models.append(model)

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

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Input observation')

        # Write header
        self._log_header1(gammalib.TERSE, 'Generate background models')

        # Clear models
        self._models.clear()

        # Generate background models
        if self['runwise'].boolean():
            self._generate_runwise_bkg()
        else:
            model = self._generate_bkg(self.obs())
            self._models.append(model)

        # Write models
        self._log_models(gammalib.NORMAL, self._models, 'Output models')

        # Return
        return

    def save(self):
        """
        Save models
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save models')

        # Get models filename
        outmodel = self['outmodel'].filename()

        # Log file name
        self._log_value(gammalib.NORMAL, 'Models file', outmodel.url())

        # Save models
        self._models.save(outmodel)

        # Return
        return

    def models(self):
        """
        Return background models

        Returns
        -------
        ebounds : `~gammalib.GModels()`
            Background models
        """
        # Return
        return self._models


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csbkgmodel(sys.argv)

    # Execute application
    app.execute()
