#! /usr/bin/env python
# ==========================================================================
# Generates the TS distribution for a particular model
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
from cscripts import modutils
from cscripts import ioutils


# ============== #
# cstsdist class #
# ============== #
class cstsdist(ctools.cslikelihood):
    """
    Generates Test Statistic distribution for a model
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Initialise application by calling the appropriate class constructor
        self._init_cslikelihood('cstsdist', ctools.__version__, argv)

        # Initialise some members
        self._srcname     = ''
        self._log_clients = False

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Set observation if not done before
        if self.obs().size() == 0:
            self.obs(self._get_observations())

        # Set observation statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

        # Get source name
        self._srcname = self['srcname'].string()

        # Query parameters for binned if requested
        if self['enumbins'].integer() != 0:
            self['npix'].integer()
            self['binsz'].real()
            self['coordsys'].string()
            self['proj'].string()

        # Set models if we have none
        if self.obs().models().size() == 0:
            self.obs().models(self['inmodel'].filename())

        # Query parameters
        self['edisp'].boolean()
        self['ntrials'].integer()
        self['outfile'].filename()
        self['debug'].boolean()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _trial(self, seed):
        """
        Create the TS for a single trial

        Parameters
        ----------
        seed : int
            Random number generator seed

        Returns
        -------
        result : dict
            Result dictionary
        """
        # Write header
        self._log_header2(gammalib.EXPLICIT, 'Trial %d' % (seed+1))

        # Set default binned parameters
        coordsys = 'CEL'
        proj     = 'TAN'
        npix     = 0
        binsz    = 0.0

        # If binned analysis is requested the read the binned parameters
        if self['enumbins'].integer() != 0:
            npix     = self['npix'].integer()
            binsz    = self['binsz'].real()
            coordsys = self['coordsys'].string()
            proj     = self['proj'].string()

        # Simulate events
        sim = obsutils.sim(self.obs(),
                           nbins = self['enumbins'].integer(),
                           seed  = seed,
                           proj  = proj,
                           coord = coordsys,
                           binsz = binsz,
                           npix  = npix,
                           log   = self._log_clients,
                           debug = self['debug'].boolean())

        # Determine number of events in simulation
        nevents = 0.0
        for run in sim:
            nevents += run.events().number()

        # Write simulation results
        self._log_header3(gammalib.EXPLICIT, 'Simulation')
        self._log_value(gammalib.EXPLICIT, 'Number of simulated events', nevents)

        # Fit model
        fit = ctools.ctlike(sim)
        fit['debug']   = self['debug'].boolean()
        fit['chatter'] = self['chatter'].integer()
        fit.run()

        # Get model fitting results
        logL   = fit.opt().value()
        npred  = fit.obs().npred()
        models = fit.obs().models()
        model  = models[self._srcname]
        ts     = model.ts()

        # Write fit results, either explicit or normal
        if self._logExplicit():
            self._log_header3(gammalib.EXPLICIT, 'Test source model fit')
            self._log_value(gammalib.EXPLICIT, 'Test statistics', ts)
            self._log_value(gammalib.EXPLICIT, 'log likelihood', logL)
            self._log_value(gammalib.EXPLICIT, 'Number of predicted events', npred)
            for model in models:
                self._log_value(gammalib.EXPLICIT, 'Model', model.name())
                for par in model:
                    self._log_string(gammalib.EXPLICIT, str(par)+'\n')
        elif self._logNormal():
            name  = 'Trial %d' % seed
            value = 'TS=%.3f  Prefactor=%e +/- %e' % \
                    (ts, model['Prefactor'].value(), model['Prefactor'].error())
            self._log_value(gammalib.TERSE, name, value)

        # Initialise results
        colnames = []
        values   = {}

        # Set TS value
        colnames.append('TS')
        values['TS'] = ts

        # Set Nevents
        colnames.append('Nevents')
        values['Nevents'] = nevents

        # Set Npred
        colnames.append('Npred')
        values['Npred'] = npred

        # Gather free full fit parameters
        for model in models:
            model_name = model.name()
            for par in model:
                if par.is_free():

                    # Set parameter name
                    name = model_name+"_"+par.name()

                    # Append value
                    colnames.append(name)
                    values[name] = par.value()

                    # Append error
                    name = 'e_'+name
                    colnames.append(name)
                    values[name] = par.error()

        # Bundle together results
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

        # Set test source model for this observation
        self.models(modutils.test_source(self.obs().models(), self._srcname))

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Input observation')

        # Write models into logger
        self._log_models(gammalib.TERSE, self.obs().models(), 'Input model')

        # Write header
        self._log_header1(gammalib.TERSE, 'Generate TS distribution')

        # Loop over trials
        for seed in range(self['ntrials'].integer()):

            # Make a trial
            result = self._trial(seed)

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
            Model container
        """
        # Copy models
        self.obs().models(models)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cstsdist(sys.argv)

    # Execute application
    app.execute()
