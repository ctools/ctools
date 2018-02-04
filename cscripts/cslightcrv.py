#! /usr/bin/env python
# ==========================================================================
# Generates a lightcurve.
#
# Copyright (C) 2014-2018 Michael Mayer
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


# ================ #
# cslightcrv class #
# ================ #
class cslightcrv(ctools.csobservation):
    """
    Generates a lightcurve

    The cslightcrv class generates a light curve for Imaging Air Cherenkov
    Telescope event data by performing a maximum likelihood fit using
    ctlike in a series of time bins. The time bins can be either
    specified in an ASCII file, as an interval divided into equally
    sized time bins, or can be taken from the Good Time Intervals of the
    observation(s).

    The format of the ASCII file is one row per time bin, each specifying
    the start of stop value of the bin, separated by a whitespace. The
    times are given in Modified Julian Days (MJD). 

    Examples:
            >>> lcrv = cslightcrv()
            >>> lcrv.run()
            >>> ... (querying for parameters) ...
            >>> fits = lcrv.lightcurve()
                Generates a light curve and retrieves the results in
                a FITS file.

            >>> lcrv = cslightcrv()
            >>> lcrv.execute()
            >>> ... (querying for parameters) ...
                Generates a light curve and saves results in a FITS file.

            >>> lcrv = cslightcrv(obs)
            >>> lcrv.execute()
            >>> ... (querying for parameters) ...
                Generates a light curve from the observations in an
                observation container and saves results in a FITS file.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise some members
        self._srcname = ''
        self._tbins   = gammalib.GGti()
        self._onoff   = False
        self._stacked = False
        self._fits    = gammalib.GFits()

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Setup observations (require response and allow event list, don't
        # allow counts cube)
        self._setup_observations(self.obs(), True, True, False)

        # Set observation statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

        # Set models if there are none in the container
        if self.obs().models().size() == 0:
            self.obs().models(self['inmodel'].filename())

        # Get source name
        self._srcname = self['srcname'].string()

        # Get time boundaries
        self._tbins = self._create_tbounds()

        # Set On/Off analysis flag and query relevant user parameters
        self._onoff = self._is_onoff()

        # Set stacked analysis flag and query relevant user parameters
        if not self._onoff:
            self._stacked = self._is_stacked()

        # Query the hidden parameters, just in case
        self['edisp'].boolean()
        self['calc_ulim'].boolean()
        self['calc_ts'].boolean()
        self['fix_bkg'].boolean()
        self['fix_srcs'].boolean()

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _create_tbounds(self):
        """
        Creates light curve time bins

        Returns
        -------
        gti : `~gammalib.GGti`
            Light curve bins in form of Good Time Intervals
        """
        # Initialise Good Time Intervals
        gti = gammalib.GGti()

        # Get algorithm to use for defining the time intervals. Possible
        # values are "FILE", "LIN" and "GTI". This is enforced at
        # parameter file level, hence no checking is needed.
        algorithm = self['tbinalg'].string()

        # If the algorithm is "FILE" then handle a FITS or an ASCII file for
        # the time bin definition
        if algorithm == 'FILE':

            # Get the filename
            filename = self['tbinfile'].filename()

            # If the file a FITS file then load GTIs directly
            if filename.is_fits():
                gti.load(filename)

            # ... otherwise load file the ASCII file as CSV file and construct
            # the GTIs from the rows of the CSV file. It is expected that the
            # CSV file has two columns containing the "START" and "STOP"
            # values of the time bins. No header row is expected.
            csv = gammalib.GCsv(filename)
            for i in range(csv.nrows()):
                tmin = gammalib.GTime()
                tmax = gammalib.GTime()
                tmin.mjd(csv.real(i,0))
                tmax.mjd(csv.real(i,1))
                gti.append(tmin,tmax)

        # If the algorithm is "LIN" then use a linear time binning, defined by
        # the "tmin", "tmax" and "tbins" user parameters
        elif algorithm == 'LIN':

            # Get start and stop time and number of time bins
            time_min = self['tmin'].time(ctools.time_reference)
            time_max = self['tmax'].time(ctools.time_reference)
            nbins    = self['tbins'].integer()

            # Compute time step in seconds and setup the GTIs
            time_step = (time_max - time_min) / float(nbins)
            for i in range(nbins):
                tmin = time_min +    i *time_step
                tmax = time_min + (i+1)*time_step
                gti.append(tmin,tmax)

        # If the algorithm is "GTI" then extract the GTIs from the observations
        # in the container and use them for the light curve time binning
        elif algorithm == 'GTI':

            # Append the GTIs of all observations
            for obs in self.obs():
                for i in range(obs.events().gti().size()):
                    gti.append(obs.events().gti().tstart(i),
                               obs.events().gti().tstop(i))

        # Return Good Time Intervals
        return gti

    def _get_free_par_names(self):
        """
        Return list of free parameter names

        Returns
        -------
        names : list of str
            List of free parameter names.
        """
        # Initialise empty list of free parameter names
        names = []

        # Collect list of free parameter names
        for par in self.obs().models()[self._srcname]:
            if par.is_free():
                names.append(par.name())

        # Return names
        return names

    def _adjust_model_pars(self):
        """
        Adjust model parameters dependent on user parameters
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Adjust model parameters')

        # Adjust model parameters dependent on input user parameters
        for model in self.obs().models():

            # Set TS flag for all models to false. The source of interest will
            # be set to true later
            model.tscalc(False)

            # Get the model classname to distinguish later sky from background
            # models
            classname = model.classname()

            # Log model name
            self._log_header3(gammalib.NORMAL, model.name())

            # If the model is the source of interest and the 'calc_ts' parameter
            # is true then enable the TS computation for the source
            if model.name() == self._srcname:
                if self['calc_ts'].boolean():
                    model.tscalc(True)

            # ... otherwise, if the model is not a sky model and the 'fix_bkg'
            # parameter is true or if the model is a sky model and the
            # 'fix_srcs' parameter is true then fix all parameters of the model
            elif ((self['fix_bkg'].boolean()  and classname != 'GModelSky') or
                  (self['fix_srcs'].boolean() and classname == 'GModelSky')):
                for par in model:
                    if par.is_free():
                        par.fix()
                        self._log_value(gammalib.NORMAL, par.name(), 'fixed')

        # Return
        return

    def _create_fits_table(self, results):
        """
        Creates FITS binary table containing light curve results

        Parameters
        ----------
        results : list of dict
            List of result dictionaries

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            FITS binary table containing light curve
        """
        # Determine number of rows in FITS table
        nrows = len(results)

        # Create FITS Table with extension "LIGHTCURVE"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('LIGHTCURVE')

        # Append time columns
        ioutils.append_result_column(table, results, 'MJD',   'days', 'mjd')
        ioutils.append_result_column(table, results, 'e_MJD', 'days', 'e_mjd')

        # Append parameter columns
        ioutils.append_model_par_column(table,
                                        self.obs().models()[self._srcname],
                                        results)

        # Append Test Statistic column "TS"
        ioutils.append_result_column(table, results, 'TS', '', 'ts')

        # Append upper limit columns
        ioutils.append_result_column(table, results, 'DiffUpperLimit',
                                     'ph/cm2/s/MeV', 'ul_diff')
        ioutils.append_result_column(table, results, 'FluxUpperLimit',
                                     'ph/cm2/s', 'ul_flux')
        ioutils.append_result_column(table, results, 'EFluxUpperLimit',
                                     'erg/cm2/s', 'ul_eflux')

        # Return table
        return table

    def _compute_ulimit(self, obs):
        """
        Computes upper limit

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Observation container

        Returns
        -------
        ul_diff, ul_flux, ul_eflux : tuple of float
            Upper limits, -1.0 of not computed
        """
        # Initialise upper flux limit
        ul_diff  = -1.0
        ul_flux  = -1.0
        ul_eflux = -1.0

        # Perform computation only if requested
        if self['calc_ulim'].boolean():

            # Write header in logger
            self._log_header3(gammalib.EXPLICIT, 'Computing upper limit')

            # Create copy of observations
            cpy_obs = obs.copy()

            # Fix parameters of source of interest in copy of observations.
            # This assures that the original spectral parameters and position
            # are used in the upper limit search. The ctulimit tool makes sure
            # that the flux parameter is freed when needed.
            source = cpy_obs.models()[self._srcname]
            for par in source:
                if par.is_free():
                    par.fix()

            # Create instance of upper limit tool
            ulimit = ctools.ctulimit(cpy_obs)
            ulimit['srcname']   = self._srcname
            ulimit['eref']      = 1.0
            ulimit['emin']      = self['emin'].real()
            ulimit['emax']      = self['emax'].real()
            ulimit['sigma_min'] = 0.0
            ulimit['sigma_max'] = 3.0

            # Try to run upper limit tool and catch any exceptions
            try:
                ulimit.run()
                ul_diff  = ulimit.diff_ulimit()
                ul_flux  = ulimit.flux_ulimit()
                ul_eflux = ulimit.eflux_ulimit()
            except:
                self._log_string(gammalib.TERSE, 'Upper limit calculation failed.')
                ul_diff  = -1.0
                ul_flux  = -1.0
                ul_eflux = -1.0

        # Return upper limit tuple
        return ul_diff, ul_flux, ul_eflux

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
        self._log_observations(gammalib.NORMAL, self.obs(), 'Observation')

        # Get time boundaries
        tmin = self._tbins.tstart(0)
        tmax = self._tbins.tstop(self._tbins.size()-1)

        # Select events
        select = ctools.ctselect(self.obs())
        select['emin'] = self['emin'].real()
        select['emax'] = self['emax'].real()
        select['tmin'] = tmin.convert(ctools.time_reference)
        select['tmax'] = tmax.convert(ctools.time_reference)
        select['rad']  = 'UNDEFINED'
        select['ra']   = 'UNDEFINED'
        select['dec']  = 'UNDEFINED'
        select.run()

        # Extract observations
        self.obs(select.obs().copy())

        # Write observation into logger
        self._log_header1(gammalib.TERSE,
                          gammalib.number('Selected observation',
                                          len(self.obs())))

        # Adjust model parameters dependent on user parameters
        self._adjust_model_pars()

        # Write header
        self._log_header1(gammalib.TERSE, 'Generate lightcurve')

        # Initialise list of result dictionaries
        results = []

        # Get source parameters
        pars = self._get_free_par_names()

        # Loop over time bins
        for i in range(self._tbins.size()):

            # Get time boundaries
            tmin = self._tbins.tstart(i)
            tmax = self._tbins.tstop(i)

            # Write time bin into header
            self._log_header2(gammalib.TERSE, 'MJD %f - %f ' %
                              (tmin.mjd(), tmax.mjd()))

            # Compute time bin center and time width
            twidth = 0.5 * (tmax - tmin) # in seconds
            tmean  = tmin + twidth

            # Initialise result dictionary
            result = {'mjd': tmean.mjd(),
                      'e_mjd': twidth / gammalib.sec_in_day,
                      'ts': 0.0,
                      'ul_diff': 0.0,
                      'ul_flux': 0.0,
                      'ul_eflux': 0.0,
                      'pars': pars,
                      'values': {}}

            # Log information
            self._log_header3(gammalib.EXPLICIT, 'Selecting events')

            # Select events
            select = ctools.ctselect(self.obs())
            select['emin'] = self['emin'].real()
            select['emax'] = self['emax'].real()
            select['tmin'] = tmin.convert(ctools.time_reference)
            select['tmax'] = tmax.convert(ctools.time_reference)
            select['rad']  = 'UNDEFINED'
            select['ra']   = 'UNDEFINED'
            select['dec']  = 'UNDEFINED'
            select.run()

            # Retrieve observation
            obs = select.obs()

            # Deal with stacked and On/Off Observations
            if self._stacked or self._onoff:

                # If a stacked analysis is requested bin the events
                # and compute the stacked response functions and setup
                # an observation container with a single stacked observation.
                if self._stacked:
                    new_obs = obsutils.get_stacked_obs(self, obs)

                # ... otherwise if On/Off analysis is requested generate
                # the On/Off observations and response
                elif self._onoff:
                    new_obs = obsutils.get_onoff_obs(self, obs)

                # Extract models
                models = new_obs.models()

                # Fix background models if required
                if self['fix_bkg'].boolean():
                    for model in models:
                        if model.classname() != 'GModelSky':
                            for par in model:
                                par.fix()

                # Put back models
                new_obs.models(models)

                # Continue with new oberservation container
                obs = new_obs

            # Header
            self._log_header3(gammalib.EXPLICIT, 'Fitting the data')

            # Do maximum likelihood model fitting
            if obs.size() > 0:
                like = ctools.ctlike(obs)
                like['edisp'] = self['edisp'].boolean()
                like.run()

                # Skip bin if no event was present
                if like.obs().logL() == 0.0:

                    # Signal skipping of bin
                    self._log_value(gammalib.TERSE, 'Warning',
                                    'No event in this time bin, skip bin.')

                    # Set all results to 0
                    for par in pars:
                        result['values'][par]      = 0.0
                        result['values']['e_'+par] = 0.0

                    # Append result
                    results.append(result)

                    # Continue with next time bin
                    continue

                # Retrieve model fitting results for source of interest
                source = like.obs().models()[self._srcname]

                # Extract parameter values
                for par in pars:
                    result['values'][par]      = source[par].value()
                    result['values']['e_'+par] = source[par].error()

                # Calculate upper limit (-1 if not computed)
                #ul_diff, ul_flux, ul_eflux = self._compute_ulimit(like.obs())
                ul_diff, ul_flux, ul_eflux = self._compute_ulimit(obs)
                if ul_diff > 0.0:
                    result['ul_diff']  = ul_diff
                    result['ul_flux']  = ul_flux
                    result['ul_eflux'] = ul_eflux

                # Extract Test Statistic value
                if self['calc_ts'].boolean():
                    result['ts'] = source.ts()

                # Append result to list of dictionaries
                results.append(result)

                # Log results for this time bin
                self._log.header3('Results')
                pars = self._get_free_par_names()
                for par in pars:
                    value = source[par].value()
                    error = source[par].error()
                    unit  = source[par].unit()
                    self._log_value(gammalib.NORMAL, par,
                                    str(value)+' +/- '+str(error)+' '+unit)
                if ul_diff > 0.0:
                    self._log_value(gammalib.NORMAL, 'Upper flux limit',
                                    str(result['ul_diff'])+' ph/cm2/s/MeV')
                    self._log_value(gammalib.NORMAL, 'Upper flux limit',
                                    str(result['ul_flux'])+' ph/cm2/s')
                    self._log_value(gammalib.NORMAL, 'Upper flux limit',
                                    str(result['ul_eflux'])+' erg/cm2/s')
                if self['calc_ts'].boolean():
                    self._log_value(gammalib.NORMAL, 'Test Statistic', result['ts'])

            # Otherwise, if observations size is 0, signal bin is skipped and
            # fill results table with zeros
            else:
                self._log_value(gammalib.TERSE, 'Warning',
                                'No observations available in this time bin, '
                                'skip bin.')

                # Set all results to 0
                for par in pars:
                    result['values'][par]        = 0.0
                    result['values']['e_' + par] = 0.0

                # Append result
                results.append(result)

                # Continue with next time bin
                continue

        # Create FITS table from results
        table = self._create_fits_table(results)

        # Create FITS file and append FITS table to FITS file
        self._fits = gammalib.GFits()
        self._fits.append(table)

        # Optionally publish light curve
        if self['publish'].boolean():
            self.publish()

        # Return
        return

    def save(self):
        """
        Save light curve
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save light curve')

        # Get light curve filename
        outfile = self['outfile'].filename()

        # Continue only filename and residual map are valid
        if self._fits != None:

            # Log file name
            self._log_value(gammalib.NORMAL, 'Light curve file', outfile.url())

            # Save spectrum
            self._fits.saveto(outfile, self._clobber())

        # Return
        return

    def publish(self, name=''):
        """
        Publish light curve

        Parameters
        ----------
        name : str, optional
            Name of light curve
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Publish light curve')

        # Continue only if light curve is valid
        if self._fits.contains('LIGHTCURVE'):

            # Set default name is user name is empty
            if not name:
                user_name = self._name()
            else:
                user_name = name

            # Log file name
            self._log_value(gammalib.NORMAL, 'Light curve name', user_name)

            # Publish light curve
            self._fits.publish('LIGHTCURVE', user_name)

        # Return
        return

    def lightcurve(self):
        """
        Return light curve FITS file

        Returns
        -------
        fits : `~gammalib.GFits()`
            FITS file containing light curve
        """
        # Return
        return self._fits

    def models(self, models):
        """
        Set model

        Parameters
        ----------
        models : `~gammalib.GModels()`
            Set model container
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
    app = cslightcrv(sys.argv)

    # Execute application
    app.execute()
