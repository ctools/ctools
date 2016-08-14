#! /usr/bin/env python
# ==========================================================================
# This script generates the pull distribution for all model parameters.
#
# Copyright (C) 2011-2016 Juergen Knoedlseder
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
from cscripts import obsutils
from cscripts import ioutils


# ============ #
# cspull class #
# ============ #
class cspull(ctools.cscript):
    """
    Generates pull distributions for a model
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'cspull'
        self._version = '1.2.0'

        # Initialise some members
        self._edisp       = False
        self._coordsys    = 'CEL'
        self._proj        = 'TAN'
        self._log_clients = False  # Static parameter
        self._offset      = 0.0
        self._ntrials     = 0
        self._npix        = 0
        self._binsz       = 0.0
        self._enumbins    = 0
        self._seed        = 1
        self._chatter     = 2

        # Initialise observation container from constructor arguments
        self._obs, argv = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # If there are no observations in container then get some ...
        if self._obs.size() == 0:
            self._obs = self._get_observations()

            # Check for requested pattern and set wobble pattern
            if self['pattern'].string() == 'four':
                self._obs = obsutils.set_observations(self['ra'].real(), self['dec'].real(), self['rad'].real(),
                                                      self['tmin'].real(), self['tmax'].real(),
                                                      self['emin'].real(), self['emax'].real(),
                                                      self['irf'].string(), self['caldb'].string(),
                                                      deadc=self['deadc'].real(),
                                                      pattern=self['pattern'].string(), offset=self['offset'].real())

        # ... otherwise add response information and energy boundaries
        # in case they are missing
        else:
            self._setup_observations(self._obs)

        # Get number of energy bins
        self._enumbins = self['enumbins'].integer()

        # Read parameters for binned if requested
        if self._enumbins != 0:
            self._npix     = self['npix'].integer()
            self._binsz    = self['binsz'].real()
            self._coordsys = self['coordsys'].string()
            self._proj     = self['proj'].string()

        # Set models if we have none
        if self._obs.models().size() == 0:
            self._obs.models(self['inmodel'].filename())

        # Read other parameters    
        self._ntrials = self['ntrials'].integer()
        self._edisp   = self['edisp'].boolean()
        self._offset  = self['offset'].real()
        self._seed    = self['seed'].integer()
        self._chatter = self['chatter'].integer()

        # Query some parameters
        self['outfile'].filename()
        self['profile'].boolean()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _obs_string(self, obs):
        """
        Generate summary string for observation

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation`
            Observation

        Returns
        -------
        text : str
            Summary string
        """
        # Extract information from observation
        emin   = obs.events().ebounds().emin().TeV()
        emax   = obs.events().ebounds().emax().TeV()
        events = obs.events().number()
        binned = (obs.events().classname() == 'GCTAEventCube')
        if binned:
            mode = 'binned'
        else:
            mode = 'unbinned'

        # Compose summary string
        if events > 0:
            text = '%d events, %.3f-%.3f TeV, %s' % (events, emin, emax, mode)
        else:
            text = '%.3f-%.3f TeV, %s' % (emin, emax, mode)

        # Return summary string
        return text

    def _trial(self, seed):
        """
        Compute the pull for a single trial

        Parameters
        ----------
        seed : int
            Random number generator seed

        Returns
        -------
        result : dict
            Dictionary of results
        """
        # Write header
        if self._logNormal():
            self._log.header2('Trial '+str(seed-self._seed+1))

        # If we have a binned obeservation then specify the lower and
        # upper energy limit in TeV
        if self._enumbins > 0:
            emin = self['emin'].real()
            emax = self['emax'].real()
        else:
            emin = None
            emax = None

        # Simulate events
        obs = obsutils.sim(self._obs,
                           emin=emin,
                           emax=emax,
                           nbins=self._enumbins,
                           addbounds=True,
                           seed=seed,
                           binsz=self._binsz,
                           npix=self._npix,
                           proj=self._proj,
                           coord=self._coordsys,
                           edisp=self._edisp,
                           log=self._log_clients,
                           debug=self._logDebug(),
                           chatter=self._chatter)

        # Determine number of events in simulation
        nevents = 0.0
        for run in obs:
            nevents += run.events().number()

        # Write simulation results
        self._log_header3(gammalib.NORMAL, 'Simulation')
        for run in self._obs:
            self._log_value(gammalib.NORMAL, 'Input observation %s' % run.id(),
                            self._obs_string(run))
        for run in obs:
            self._log_value(gammalib.NORMAL, 'Output observation %s' % run.id(),
                            self._obs_string(run))
        self._log_value(gammalib.NORMAL, 'Number of simulated events', nevents)

        # Fit model
        if self['profile'].boolean():
            models = self._obs.models()
            for model in models:
                like = ctools.cterror(obs)
                like['srcname'] = model.name()
                like['edisp']   = self._edisp
                like['debug']   = self._logDebug()
                like['chatter'] = self._chatter
                like.run()
        else:
            like = ctools.ctlike(obs)
            like['edisp']   = self._edisp
            like['debug']   = self._logDebug()
            like['chatter'] = self._chatter
            like.run()

        # Store results
        logL   = like.opt().value()
        npred  = like.obs().npred()
        models = like.obs().models()

        # Write result header
        if self._logNormal():
            self._log.header3('Pulls')

        # Gather results in form of a list of result columns and a
        # dictionary containing the results. The result contains the
        # log-likelihood, the number of simulated events, the number of
        # predicted events and for each fitted parameter the fitted value,
        # the pull and the fit error.
        #
        # Note that we do not use the model and parameter iterators
        # because we need the indices to get the true (or real) parameter
        # values from the input models.
        colnames = ['LogL', 'Sim_Events', 'Npred_Events']
        values   = {'LogL': logL, 'Sim_Events': nevents, 'Npred_Events': npred}
        for i in range(models.size()):
            model = models[i]
            for k in range(model.size()):
                par = model[k]
                if par.is_free():

                    # Set name as a combination of model name and parameter
                    # name separated by an underscore. In that way each
                    # parameter has a unique name.
                    name = model.name()+'_'+par.name()

                    # Append parameter, Pull_parameter and e_parameter column
                    # names
                    colnames.append(name)
                    colnames.append('Pull_'+name)
                    colnames.append('e_'+name)

                    # Compute pull for this parameter as the difference
                    #               (fitted - true) / error
                    # In case that the error is 0 the pull is set to 99
                    fitted_value = par.value()
                    real_value   = self._obs.models()[i][k].value()
                    error        = par.error()
                    if error != 0.0:
                        pull = (fitted_value - real_value) / error
                    else:
                        pull = 99.0

                    # Store results in dictionary
                    values[name]         = fitted_value
                    values['Pull_'+name] = pull
                    values['e_'+name]    = error

                    # Write results into logger
                    value = '%.4f (%e +/- %e)' % (pull, fitted_value, error)
                    self._log_value(gammalib.NORMAL, name, value)

        # Bundle together results in a dictionary
        result = {'colnames': colnames, 'values': values}

        # Return
        return result


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
        self._log_observations(gammalib.NORMAL, self._obs, 'Input observation')

        # Write header
        self._log_header1(gammalib.TERSE, 'Generate pull distribution')

        # Loop over trials
        for seed in range(self._ntrials):

            # Make a trial and add initial seed
            result = self._trial(seed + self._seed)

            # Write out trial result
            ioutils.write_csv_row(self['outfile'].filename().url(), seed,
                                  result['colnames'], result['values'])

        # Return
        return

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Run the script
        self.run()

        # Return
        return

    def models(self, models):
        """
        Set model

        Parameters
        ----------
        models : `~gammalib.GModels`
            Set model container
        """
        # Set model container
        self._obs.models(models)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cspull(sys.argv)

    # Execute application
    app.execute()
