#! /usr/bin/env python
# ==========================================================================
# This script generates the pull distribution for all model parameters.
#
# Copyright (C) 2011-2017 Juergen Knoedlseder
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
from cscripts import obsutils
from cscripts import ioutils


# ============ #
# cspull class #
# ============ #
class cspull(ctools.csobservation):
    """
    Generates pull distributions for a model
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # If there are no observations in container then get some ...
        if self.obs().is_empty():
            self.obs(self._get_observations())

        # ... otherwise add response information and energy boundaries
        # in case they are missing
        else:
            self._setup_observations(self.obs())

        # Set observation statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

        # Get number of energy bins
        enumbins = self['enumbins'].integer()

        # Query parameters for On/Off observation
        if len(self['onsrc'].string()) > 0:
            self['onrad'].real()

        # Query parameters for binned if requested
        elif enumbins > 0:
            self['npix'].integer()
            self['binsz'].real()
            self['coordsys'].string()
            self['proj'].string()

        # Set models if we have none
        if self.obs().models().is_empty():
            self.obs().models(self['inmodel'].filename())

        # Query other parameters
        self['ntrials'].integer()
        self['edisp'].boolean()
        self['seed'].integer()
        self['chatter'].integer()

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
        if obs.classname() == 'GCTAOnOffObservation':
            emin = obs.on_spec().ebounds().emin().TeV()
            emax = obs.on_spec().ebounds().emin().TeV()
            mode = 'On/Off'
        else:
            emin = obs.events().ebounds().emin().TeV()
            emax = obs.events().ebounds().emax().TeV()
            binned = (obs.events().classname() == 'GCTAEventCube')
            if binned:
                mode = 'binned'
            else:
                mode = 'unbinned'
        events = obs.nobserved()

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
        self._log_header2(gammalib.NORMAL, 'Trial %d' %
                          (seed-self['seed'].integer()+1))

        # Get number of energy bins and On source name and initialise
        # some parameters
        nbins    = self['enumbins'].integer()
        onsrc    = self['onsrc'].string()
        edisp    = self['edisp'].boolean()
        emin     = None
        emax     = None
        binsz    = 0.0
        npix     = 0
        proj     = 'TAN'
        coordsys = 'CEL'

        # If we have a On source name then set On region radius
        if len(onsrc) > 0:
            onrad = self['onrad'].real()
            emin  = self['emin'].real()
            emax  = self['emax'].real()
            edisp = True   # Use always energy dispersion for On/Off
        else:

            # Reset On region source name and radius
            onrad = 0.0
            onsrc = None

            # If we have a binned obeservation then specify the lower and
            # upper energy limit in TeV
            if nbins > 0:
                emin     = self['emin'].real()
                emax     = self['emax'].real()
                binsz    = self['binsz'].real()
                npix     = self['npix'].integer()
                proj     = self['proj'].string()
                coordsys = self['coordsys'].string()

        # Simulate events
        obs = obsutils.sim(self.obs(),
                           emin=emin, emax=emax, nbins=nbins,
                           onsrc=onsrc, onrad=onrad,
                           addbounds=True, seed=seed,
                           binsz=binsz, npix=npix, proj=proj, coord=coordsys,
                           edisp=edisp, log=False, debug=self._logDebug(),
                           chatter=self['chatter'].integer())

        # Determine number of events in simulation
        nevents = 0.0
        for run in obs:
            nevents += run.nobserved()

        # Write simulation results
        self._log_header3(gammalib.NORMAL, 'Simulation')
        for run in self.obs():
            self._log_value(gammalib.NORMAL, 'Input observation %s' % run.id(),
                            self._obs_string(run))
        for run in obs:
            self._log_value(gammalib.NORMAL, 'Output observation %s' % run.id(),
                            self._obs_string(run))
        self._log_value(gammalib.NORMAL, 'Number of simulated events', nevents)

        # Fit model
        if self['profile'].boolean():
            models = self.obs().models()
            for model in models:
                like = ctools.cterror(obs)
                like['srcname'] = model.name()
                like['edisp']   = edisp
                like['debug']   = self._logDebug()
                like['chatter'] = self['chatter'].integer()
                like.run()
        else:
            like = ctools.ctlike(obs)
            like['edisp']   = edisp
            like['debug']   = self._logDebug()
            like['chatter'] = self['chatter'].integer()
            like.run()

        # Store results
        logL   = like.opt().value()
        npred  = like.obs().npred()
        models = like.obs().models()

        # Write result header
        self._log_header3(gammalib.NORMAL, 'Pulls')

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
                    real_value   = self.obs().models()[i][k].value()
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
        self._log_observations(gammalib.NORMAL, self.obs(), 'Input observation')

        # Write header
        self._log_header1(gammalib.TERSE, 'Generate pull distribution')

        # Loop over trials
        for seed in range(self['ntrials'].integer()):

            # Make a trial and add initial seed
            result = self._trial(seed + self['seed'].integer())

            # Write out trial result
            ioutils.write_csv_row(self['outfile'].filename().url(), seed,
                                  result['colnames'], result['values'])

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
        self.obs().models(models)

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
