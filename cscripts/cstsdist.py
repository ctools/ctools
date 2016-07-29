#! /usr/bin/env python
# ==========================================================================
# Generates the TS distribution for a particular model
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
import csv
import gammalib
import ctools
from cscripts import obsutils


# ============== #
# cstsdist class #
# ============== #
class cstsdist(ctools.cscript):
    """
    Generates Test Statistic distribution for a model
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = 'cstsdist'
        self._version = '1.2.0'

        # Initialise some members
        self._srcname     = ''
        self._log_clients = False

        # Initialise observation container from constructor arguments
        self._obs, argv = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Set observation if not done before
        if self._obs == None or self._obs.size() == 0:
            self._obs = self._get_observations()

            # Check for requested pattern and use above observation parameters
            # to set wobble pattern
            if self['pattern'].string() == 'four':
                self._obs = self._set_obs()

        # Get source name
        self._srcname = self['srcname'].string()

        # Query parameters for binned if requested
        if self['enumbins'].integer() != 0:
            self['npix'].integer()
            self['binsz'].real()
            self['coordsys'].string()
            self['proj'].string()

        # Set models if we have none
        if self._obs.models().size() == 0:
            self._obs.models(self['inmodel'].filename())

        # Query parameters
        self['edisp'].boolean()
        self['ntrials'].integer()
        self['outfile'].filename()
        self['debug'].boolean()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log('\n')

        # Return
        return

    def _set_obs(self):
        """
        Set an observation container

        Returns
        -------
        obs : `~gammalib.GObservations`
            Observation container
        """
        # Set duration
        duration = self['tmax'].real() - self['tmin'].real()

        # Setup observation definition list
        obsdeflist = obsutils.set_obs_patterns(self['pattern'].string(),
                                               ra     = self['ra'].real(),
                                               dec    = self['dec'].real(),
                                               offset = self['offset'].real())

        # Create list of observations
        obs = obsutils.set_obs_list(obsdeflist,
                                    tstart   = self['tmin'].real(),
                                    duration = duration,
                                    deadc    = self['deadc'].real(),
                                    emin     = self['emin'].real(),
                                    emax     = self['emax'].real(),
                                    rad      = self['rad'].real(),
                                    irf      = self['irf'].string(),
                                    caldb    = self['caldb'].string())

        # Return observation container
        return obs

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
        if self._logExplicit():
            self._log.header2('Trial '+str(seed+1))

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
        sim = obsutils.sim(self._obs,
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
        if self._logExplicit():
            self._log.header3('Simulation')
            self._log.parformat('Number of simulated events')
            self._log(nevents)
            self._log('\n')

        # Fit model
        fit = obsutils.fit(sim, log   = self._log_clients,
                                debug = self['debug'].boolean())

        # Get model fitting results
        logL   = fit.opt().value()
        npred  = fit.obs().npred()
        models = fit.obs().models()
        model  = models[self._srcname]
        ts     = model.ts()

        # Write fit results, either explicit or normal
        if self._logExplicit():
            self._log.header3('Test source model fit')
            self._log.parformat('Test statistics')
            self._log(ts)
            self._log('\n')
            self._log.parformat('log likelihood')
            self._log(logL)
            self._log('\n')
            self._log.parformat('Number of predicted events')
            self._log(npred)
            self._log('\n')
            for model in models:
                self._log.parformat('Model')
                self._log(model.name())
                self._log('\n')
                for par in model:
                    self._log(str(par)+'\n')
        elif self._logNormal():
            self._log.parformat('Trial '+str(seed))
            self._log('TS=')
            self._log(ts)
            self._log('  Prefactor=')
            self._log(model['Prefactor'].value())
            self._log('+/-')
            self._log(model['Prefactor'].error())
            self._log('\n')

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

        # Set models for this observation
        self.models(obsutils.set_ts_model(self._obs.models(), self._srcname))

        # Write models into logger
        if self._logTerse():
            self._log('\n')
            self._log.header1('Models')
            self._log(str(self._obs.models()))
            self._log('\n')

        # Write observation into logger
        if self._logTerse():
            self._log('\n')
            self._log.header1(gammalib.number('Observation',len(self._obs)))
            self._log(str(self._obs))
            self._log('\n')
        if self._logExplicit():
            for obs in self._obs:
                self._log(str(obs))
                self._log('\n')

        # Write header
        if self._logTerse():
            self._log('\n')
            self._log.header1('Generate TS distribution')

        # Loop over trials
        for seed in range(self['ntrials'].integer()):

            # Make a trial
            result = self._trial(seed)

            # Set output filename
            outfile = self['outfile'].filename().url()

            # Write out result immediately
            if seed == 0:
                f      = open(outfile, 'w')
                writer = csv.DictWriter(f, result['colnames'])
                headers = {}
                for n in result['colnames']:
                    headers[n] = n
                writer.writerow(headers)
            else:
                f = open(outfile, 'a')
            writer = csv.DictWriter(f, result['colnames'])
            writer.writerow(result['values'])
            f.close()

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
        self._obs.models(models)

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


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cstsdist(sys.argv)

    # Execute application
    app.execute()
