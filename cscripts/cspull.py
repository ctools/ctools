#! /usr/bin/env python
# ==========================================================================
# This script generates the pull distribution for all model parameters.
#
# Copyright (C) 2011-2021 Juergen Knoedlseder
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
from cscripts import mputils


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

        # Initialise class members
        self._fits     = None
        self._nthreads = 0

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
        state = {'base'     : ctools.csobservation.__getstate__(self),
                 'fits'     : self._fits,
                 'nthreads' : self._nthreads}

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
        self._fits     = state['fits']
        self._nthreads = state['nthreads']

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
        if gammalib.toupper(self['onsrc'].string()) != 'NONE':
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
        self['profile'].boolean()

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Set number of processes for multiprocessing
        self._nthreads = mputils.nthreads(self)

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
            emax = obs.on_spec().ebounds().emax().TeV()
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
        nbins     = self['enumbins'].integer()
        onsrc     = self['onsrc'].string()
        edisp     = self['edisp'].boolean()
        statistic = self['statistic'].string()
        emin      = None
        emax      = None
        binsz     = 0.0
        npix      = 0
        proj      = 'TAN'
        coordsys  = 'CEL'

        # If we have a On source name then set On region radius
        if gammalib.toupper(onsrc) != 'NONE':
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
                           edisp=edisp, nthreads=1, log=False,
                           debug=self._logDebug(), chatter=self['chatter'].integer())

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
                like['srcname']   = model.name()
                like['edisp']     = edisp
                like['statistic'] = statistic
                like['nthreads']  = 1  # Avoids OpenMP conflict
                like['debug']     = self._logDebug()
                like['chatter']   = self['chatter'].integer()
                like.run()
        else:
            like = ctools.ctlike(obs)
            like['edisp']     = edisp
            like['statistic'] = statistic
            like['nthreads']  = 1  # Avoids OpenMP conflict
            like['debug']     = self._logDebug()
            like['chatter']   = self['chatter'].integer()
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
                    name = model.name() + '_' + par.name()

                    # Set "Parameter", "Parameter_error" and "Parameter_Pull"
                    # column names
                    col_par       = name
                    col_par_error = name+'_error'
                    col_par_pull  = name+'_pull'

                    # Append column names
                    colnames.append(col_par)
                    colnames.append(col_par_error)
                    colnames.append(col_par_pull)

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
                    values[col_par]       = fitted_value
                    values[col_par_error] = error
                    values[col_par_pull]  = pull

                    # Write results into logger
                    value = '%.4f (%e +/- %e)' % (pull, fitted_value, error)
                    self._log_value(gammalib.NORMAL, name, value)

        # Bundle together results in a dictionary
        result = {'colnames': colnames, 'values': values}

        # Return
        return result

    def _create_fits(self, results):
        """
        Create FITS file from results

        Parameters
        ----------
        results : list of dict
            List of result dictionaries
        """
        # Gather headers for parameter columns
        headers = []
        for colname in results[0]['colnames']:
            if colname != 'LogL' and colname != 'Sim_Events' and \
               colname != 'Npred_Events':
                headers.append(colname)

        # Create FITS table columns
        nrows   = len(results)
        logl    = gammalib.GFitsTableDoubleCol('LOGL', nrows)
        nevents = gammalib.GFitsTableDoubleCol('NEVENTS', nrows)
        npred   = gammalib.GFitsTableDoubleCol('NPRED', nrows)
        logl.unit('')
        nevents.unit('counts')
        npred.unit('counts')
        columns = []
        for header in headers:
            name   = gammalib.toupper(header)
            column = gammalib.GFitsTableDoubleCol(name, nrows)
            column.unit('')
            columns.append(column)

        # Fill FITS table columns
        for i, result in enumerate(results):
            logl[i]    = result['values']['LogL']
            nevents[i] = result['values']['Sim_Events']
            npred[i]   = result['values']['Npred_Events']
            for k, column in enumerate(columns):
                column[i] = result['values'][headers[k]]

        # Initialise FITS Table with extension "PULL_DISTRIBUTION"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('PULL_DISTRIBUTION')

        # Add keywords for compatibility with gammalib.GMWLSpectrum
        table.card('INSTRUME', 'CTA', 'Name of Instrument')
        table.card('TELESCOP', 'CTA', 'Name of Telescope')

        # Stamp header
        self._stamp(table)

        # Set optional keyword values
        if gammalib.toupper(self['onsrc'].string()) != 'NONE':
            onrad = self['onrad'].real()
        else:
            onrad = 0.0
        if self['enumbins'].integer() > 0:
            npix     = self['npix'].integer()
            binsz    = self['binsz'].real()
            coordsys = self['coordsys'].string()
            proj     = self['proj'].string()
        else:
            npix     = 0
            binsz    = 0.0
            coordsys = ''
            proj     = ''

        # Add script keywords
        table.card('NPULLS',   self['ntrials'].integer(),  'Number of pulls')
        table.card('SEED',     self['seed'].integer(),     'Seed value for pulls')
        table.card('ONSRC',    self['onsrc'].string(),     'Name of On surce for On/Off analysis')
        table.card('ONRAD',    onrad,                      '[deg] Radius of On region')
        table.card('ENUMBINS', self['enumbins'].integer(), 'Number of energy bins')
        table.card('NPIX',     npix,                       'Number of pixels for binned analysis')
        table.card('BINSZ',    binsz,                      'Pixel size for binned analysis')
        table.card('COORDSYS', coordsys,                   'Coordinate system for binned analysis')
        table.card('PROJ',     proj,                       'Projection for binned analysis')
        table.card('STAT',     self['statistic'].string(), 'Optimization statistic')
        table.card('EDISP',    self['edisp'].boolean(),    'Use energy dispersion?')
        table.card('PROFILE',  self['profile'].boolean(),  'Use likelihood profile method for errors?')

        # Append filled columns to fits table
        table.append(logl)
        table.append(nevents)
        table.append(npred)
        for column in columns:
            table.append(column)

        # Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append(table)

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
        self._log_header1(gammalib.TERSE, 'Generate pull distribution')

        # Get number of trials
        ntrials = self['ntrials'].integer()

        # Get seed value
        seed = self['seed'].integer()

        # Initialise results
        results  = []

        # If more than a single thread is requested then use multiprocessing
        if self._nthreads > 1:
            args        = [(self, '_trial', i + seed) for i in range(ntrials)]
            poolresults = mputils.process(self._nthreads, mputils.mpfunc, args)

        # Continue with regular processing
        for i in range(ntrials):

            # If multiprocessing was used then recover results and put them
            # into the log file
            if self._nthreads > 1:
                results.append(poolresults[i][0])
                self._log_string(gammalib.TERSE, poolresults[i][1]['log'], False)

            # ... otherwise make a trial now
            else:

                # Run trial
                result = self._trial(i + seed)

                # Append results
                results.append(result)

        # Create FITS file
        self._create_fits(results)

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

    def save(self):
        """
        Save pull FITS file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save pull distribution')

        # Continue only if FITS file is valid
        if self._fits != None:

            # Get outmap parameter
            outfile = self['outfile'].filename()

            # Log file name
            self._log_value(gammalib.NORMAL, 'Pull distribution file', outfile.url())

            # Save pull distribution
            self._fits.saveto(outfile, self['clobber'].boolean())

        # Return
        return

    def pull_distribution(self):
        """
        Return pull distribution FITS file

        Returns:
            FITS file containing pull distribution
        """
        # Return
        return self._fits


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cspull(sys.argv)

    # Execute application
    app.execute()
