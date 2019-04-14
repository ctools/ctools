#! /usr/bin/env python
# ==========================================================================
# Generates background model
#
# Copyright (C) 2018-2019 Juergen Knoedlseder
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
        state = {'base'       : ctools.csobservation.__getstate__(self),
                 'models'     : self._models,
                 'instrument' : self._instrument}

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
        self._models     = state['models']
        self._instrument = state['instrument']

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Setup observations (require response, allow event list and counts
        # cube)
        self._setup_observations(self.obs(), True, True, True)

        # Set instrument name
        self._instrument = self._get_instrument()

        # Query input parameters
        spatial = self['spatial'].string()
        if spatial == 'GAUSS(E)':
            snumbins = self['snumbins'].integer()
            if snumbins > 1:
                self['smin'].real()
                self['smax'].real()
        if spatial == 'GAUSS' or spatial == 'GAUSS(E)':
            self['gradient'].boolean()
        spectral = self['spectral'].string()
        if spectral == 'NODES':
            self._create_energies()
        self['runwise'].boolean()
        self['emin'].real()
        self['emax'].real()
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

    def _generate_initial_model(self, sigma_min=2.0, sigma_max=1000.0):
        """
        Generate initial background model

        Parameters
        ----------
        sigma_min : float, optional
            Minimum sigma
        sigma_max : float, optional
            Maximum sigma

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

            # Set spatial model
            factor1 = gammalib.GCTAModelRadialGauss(3.0)
            factor1['Sigma'].min(sigma_min)
            factor1['Sigma'].max(sigma_max)

            # Optionally add gradient
            if self['gradient'].boolean():
                spatial = gammalib.GCTAModelSpatialMultiplicative()
                factor2 = gammalib.GCTAModelSpatialGradient()
                spatial.append(factor1)
                spatial.append(factor2)
            else:
                spatial = factor1

            # Set spectral model
            epivot   = gammalib.GEnergy(1.0, 'TeV')
            spectral = gammalib.GModelSpectralPlaw(3.0e-4, -1.5, epivot)
            spectral['Prefactor'].min(1.0e-8)

            # Set background model
            model = gammalib.GCTAModelBackground(spatial, spectral)

        # Handle GAUSS(E) model
        elif self['spatial'].string() == 'GAUSS(E)':

            # Set spatial model
            if self['snumbins'].integer() == 1:
                factor1 = gammalib.GCTAModelRadialGauss(3.0)
                factor1['Sigma'].min(sigma_min)
                factor1['Sigma'].max(sigma_max)
            else:
                emin     = gammalib.GEnergy(self['smin'].real(), 'TeV')
                emax     = gammalib.GEnergy(self['smax'].real(), 'TeV')
                energies = gammalib.GEnergies(self['snumbins'].integer(), emin, emax)
                spectrum = gammalib.GModelSpectralConst(3.0)
                nodes    = gammalib.GModelSpectralNodes(spectrum, energies)
                for i in range(nodes.nodes()):
                    nodes[i*2+1].min(sigma_min)
                    nodes[i*2+1].max(sigma_max)
                nodes.autoscale()
                factor1 = gammalib.GCTAModelSpatialGaussSpectrum(nodes)

            # Optionally add gradient
            if self['gradient'].boolean():
                spatial = gammalib.GCTAModelSpatialMultiplicative()
                factor2 = gammalib.GCTAModelSpatialGradient()
                spatial.append(factor1)
                spatial.append(factor2)
            else:
                spatial = factor1

            # Set spectral model
            epivot   = gammalib.GEnergy(1.0, 'TeV')
            spectral = gammalib.GModelSpectralPlaw(3.0e-4, -1.5, epivot)
            spectral['Prefactor'].min(1.0e-8)

            # Set background model
            model = gammalib.GCTAModelBackground(spatial, spectral)

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
        select['emin']     = self['emin'].real()
        select['emax']     = self['emax'].real()
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

        # Define node energies
        energies  = self._create_energies()

        # Create node spectrum
        nodes = gammalib.GModelSpectralNodes(spectrum, energies)
        nodes.autoscale()

        # Set minimum node intensities (this avoids NaNs during model fitting)
        for i in range(nodes.size()):
            nodes[i].min(1.0e-10*nodes[i].value())

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

        # If a NODES model is requested then refit a node spectrum
        if self['spectral'].string() == 'NODES':

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

            # Remove nodes with zero errors as they are not constrained
            # by the data and may lead to fitting problems later
            spectral = model.spectral()
            nodes    = spectral.nodes()
            for i in range(nodes):
                iint = 2*(nodes - i) - 1
                if spectral[iint].error() == 0.0:
                    spectral.remove(nodes - i - 1)

        # Return model
        return model

    def _generate_runwise_bkg(self):
        """
        Generate background models
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
