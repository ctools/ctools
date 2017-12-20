#! /usr/bin/env python
# ==========================================================================
# Generates a residual spectrum.
#
# Copyright (C) 2017- Luigi Tibaldo
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
import gammalib
import ctools
from cscripts import obsutils
import sys


# =============== #
# csresspec class #
# =============== #
class csresspec(ctools.csobservation):
    """
    Generates a residual map
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__,
                                 argv)

        # Initialise class members
        self._use_maps = False
        self._stack = False
        self._mask = False
        self._fits = gammalib.GFits()

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Setup observations (require response and allow event list as well as
        # counts cube)
        self._setup_observations(self.obs(), True, True, True)

        # Set observation statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

        # Collect number of unbinned, binned and On/Off observations in
        # observation container
        n_unbinned = 0
        n_binned = 0
        n_onoff = 0
        for obs in self.obs():
            if obs.classname() == 'GCTAObservation':
                if obs.eventtype() == 'CountsCube':
                    n_binned += 1
                else:
                    n_unbinned += 1
            elif obs.classname() == 'GCTAOnOffObservation':
                n_onoff += 1
        n_cta = n_unbinned + n_binned + n_onoff
        n_other = self.obs().size() - n_cta

        # Log census of input observations
        self._log_value(gammalib.TERSE, 'Unbinned CTA observations', n_unbinned)
        self._log_value(gammalib.TERSE, 'Binned CTA observations', n_binned)
        self._log_value(gammalib.TERSE, 'On/off CTA observations', n_onoff)
        self._log_value(gammalib.TERSE, 'Other observations', n_other)
        if n_other > 0:
            msg = 'WARNING: Only CTA observation can be handled, all non-CTA ' \
                  + 'observations will be ignored.\n'
            self._log_string(gammalib.TERSE, msg)

        # Query wether to compute model for individual components
        components = self['components'].boolean()

        # If there is only one binned observation and no model for individual
        # components is required, query for precomputed model file and set
        # use_maps to True
        if self.obs().size() == 1 and n_binned == 1 and not components:
            modcube = self['modcube'].filename()
            if modcube != 'NONE':
                self._use_maps = True

        # If there are unbinned observations query the energy binning parameters
        if n_unbinned != 0:
            self['ebinalg'].string()
            if self['ebinalg'].string() == 'FILE':
                self['ebinfile'].filename()
            else:
                self['emin'].real()
                self['emax'].real()
                self['enumbins'].integer()
            msg = 'User defined energy binning will be used for %s unbinned observations.' % n_unbinned
            self._log_string(gammalib.TERSE, msg)
            if n_cta > n_unbinned:
                n_notunbin = n_cta - n_unbinned
                msg = 'The intrinsic binning will be used for the remaining %s CTA observations.' % (
                    n_notunbin)
                self._log_string(gammalib.TERSE, msg)

        # If there is more than one observation, and observations are all
        # unbinned or all onoff query user to know if they wish stacked results
        if self.obs().size() > 1 and (
                        n_unbinned == self.obs().size() or n_onoff == self.obs().size()):
            self._stack = self['stack'].boolean()
            # If we are to stack event lists query parameters for cube creation
            if self._stack and n_unbinned == self.obs().size():
                self['coordsys'].string()
                self['proj'].string()
                self['xref'].real()
                self['yref'].real()
                self['nxpix'].integer()
                self['nypix'].integer()
                self['binsz'].real()

        # If we are not using a precomputed model and no models are available
        # in the observation container query input XML model file
        if not self._use_maps and self.obs().models().size() == 0:
            self.obs().models(self['inmodel'].filename())

        # Unless all observations are On/Off query for mask definition
        if n_onoff == n_cta:
            pass
        else:
            self._mask = self['mask'].boolean()
            if self._mask:
                self['ra'].real()
                self['dec'].real()
                self['rad'].real()
                self['regfile'].filename()

        # Unless all observations are On/Off, or we are using precomputed model
        # maps query whether to use energy dispersion
        if n_onoff == n_cta or self._use_maps:
            msg = 'Energy dispersion is applied based on the input data/model ' \
                  + 'and not according to the edisp parameter'
            self._log_string(gammalib.TERSE, msg)
        else:
            self['edisp'].boolean()

        # Query algorithm for residual computation
        self['algorithm'].string()

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()

        # Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _stack_observations(self):
        """
        Stack multiple observations and replace observation container
        with single stacked observation

        :return: `~gammalib.GObservations' container with a single stacked observation
        """
        # Use first observation to determine the type and apply the
        # appropriate stacking method

        # If On/Off use the GCTAOnOffObservation constructor
        if self.obs()[0].classname() == 'GCTAOnOffObservation':
            msg = 'Stacking %s On/Off observations.' % (self.obs().size())
            self._log_string(gammalib.TERSE, msg)
            stacked_obs = gammalib.GCTAOnOffObservation(self.obs())

        # If event list then bin observations
        elif self.obs()[0].classname() == 'GCTAObservation':
            msg = 'Stacking %s event lists.' % (self.obs().size())
            self._log_string(gammalib.TERSE, msg)
            stacked_obs = obsutils.get_stacked_obs(self, self.obs())

        # Return
        return stacked_obs

    def _bin_evlist(self, obs):
        """
        Turn single event list into counts cube
        :param   obs: `~gammalib.GObservations' observation container with single event list
        :return: `~gammalib.GObservations' binned observation container
                 dict with event list ROI and energy range information

        """
        # Retrieve information about ROI in event list
        roi = obs[0].roi()
        ra = roi.centre().dir().ra_deg()
        dec = roi.centre().dir().dec_deg()
        rad = roi.radius()

        # We will cover the whole ROI with 0.02 deg binning
        npix = int(2 * rad / 0.02) + 1

        # Bin events
        cntcube = ctools.ctbin(obs)
        cntcube['xref'] = ra
        cntcube['yref'] = dec
        cntcube['binsz'] = 0.02
        cntcube['nxpix'] = npix
        cntcube['nypix'] = npix
        cntcube['proj'] = 'TAN'
        cntcube['coordsys'] = 'CEL'
        cntcube['ebinalg'] = self['ebinalg'].string()
        if self['ebinalg'].string() == 'FILE':
            cntcube['ebinfile'] = self['ebinfile'].filename()
        else:
            cntcube['enumbins'] = self['enumbins'].integer()
            cntcube['emin'] = self['emin'].real()
            cntcube['emax'] = self['emax'].real()
        cntcube.run()

        # Retrieve the binned observation container
        binned_obs = cntcube.obs().copy()

        # Compute binned response
        response = obsutils.get_stacked_response(obs, ra, dec,
                                                 binsz=0.02,
                                                 nxpix=npix,
                                                 nypix=npix,
                                                 coordsys='CEL',
                                                 proj='TAN',
                                                 emin=self['emin'].real(),
                                                 emax=self['emax'].real(),
                                                 edisp=self['edisp'].boolean())

        if self['edisp'].boolean():
            binned_obs[0].response(response['expcube'], response['psfcube'],
                                   response['edispcube'], response['bkgcube'])
        else:
            binned_obs[0].response(response['expcube'], response['psfcube'],
                                   response['bkgcube'])

        # Get new models
        models = response['models']

        # Set models for observation container
        binned_obs.models(models)

        # Check if energy boundaries provided by user extend beyond
        # the content of the event list
        if self['emin'].real() > obs[0].events().emin().TeV():
            emin = 'INDEF'
        else:
            emin = obs[0].events().emin().TeV()
        if self['emax'].real() < obs[0].events().emax().TeV():
            emax = 'INDEF'
        else:
            emax = obs[0].events().emax().TeV()

        # Put ROI and E bound info in dictionary
        info = {'was_list': True, 'roi_ra': ra, 'roi_dec': dec, 'roi_rad': rad,
                'emin': emin, 'emax': emax}

        # Return new oberservation container
        return binned_obs, info

    def _masked_cube(self, cube, ra, dec, rad, emin='INDEF', emax='INDEF',
                     regfile='NONE'):
        """
        Mask an event cube and returns the masked cube
        :param cube: `~gammalib.GCTAEventCube' Events cube
        :param ra: double (str 'INDEF' for no selection on direction)
        :param dec: double (str 'INDEF' for no selection on direction)
        :param rad: double (str 'INDEF' for no selection on direction)
        :param emin: double (str 'INDEF' for no selection on energy)
        :param emax: double (str 'INDEF' for no selection on energy)
        :param regfile: str with name of ds9 file with exclusion regions
        :return: `~gammalib.GCTAEventCube'
        """

        # Turn cube into observation container to feed to ctcubemask
        obs = gammalib.GCTAObservation()
        obs.events(cube)
        obs_cont = gammalib.GObservations()
        obs_cont.append(obs)

        # Use ctcubemask to mask
        cubemask = ctools.ctcubemask(obs_cont)
        cubemask['ra'] = ra
        cubemask['dec'] = dec
        cubemask['rad'] = rad
        cubemask['emin'] = emin
        cubemask['emax'] = emax
        cubemask['regfile'] = regfile
        cubemask.run()

        # Return
        # Note: we return a copy to avoid memory leaks in SWIG
        return cubemask.cube().copy()

    def _cube_to_spectrum(self, cube, evlist_info):
        """
        Derive from event cube a count spectrum. If data come from event list
        use only the ROI and energy range of the original data. Apply user
        defined mask if requested.
        :param cube: `~gammalib.GCTAEventCube' input cube
        :param evlist_info: dic with information on original event list
        :return: `~gammalib.GNdarray' count spectrum
        """
        # If we started from event list mask the ROI only
        # for residual computation
        if evlist_info['was_list']:
            cube = self._masked_cube(cube, evlist_info['roi_ra'],
                                     evlist_info['roi_dec'],
                                     evlist_info['roi_rad'],
                                     emin=evlist_info['emin'],
                                     emax=evlist_info['emax'])

        # Apply user mask
        if self._mask:
            cube = self._masked_cube(cube, self['ra'], self['dec'],
                                     self['rad'],
                                     regfile=self['regfile'])

        # Extract skymap and clip at 0 to null masked areas
        counts = cube.counts().copy()
        counts = counts.clip(0.)

        # Convert skymap into GNdarray count spectrum
        counts = counts.total_counts()

        # Return
        return counts

    def _residuals_table(self, obs_id, ebounds, counts, model, residuals):
        """
        Create a Fits Table and store counts, model, and residuals
        :param obs_id: str observation id
        :param ebounds: `~gammalib.GEbounds' energy bounds
        :param counts: `~gammalib.GNdarray'
        :param model: `~gammalib.GNdarray'
        :param residuals: `~gammalib.GNdarray'
        :return: `~gammalib.GFitsBinTable'
        """
        # Create FITS table columns
        nrows = ebounds.size()
        energy_low = gammalib.GFitsTableDoubleCol('ed_Energy', nrows)
        energy_high = gammalib.GFitsTableDoubleCol('eu_Energy', nrows)
        counts_col = gammalib.GFitsTableDoubleCol('Counts', nrows)
        model_col = gammalib.GFitsTableDoubleCol('Model', nrows)
        resid_col = gammalib.GFitsTableDoubleCol('Residuals', nrows)
        energy_low.unit('TeV')
        energy_high.unit('TeV')

        # Fill FITS table columns
        for i in range(nrows):
            energy_low[i] = ebounds.emin(i).TeV()
            energy_high[i] = ebounds.emax(i).TeV()
            counts_col[i] = counts[i]
            model_col[i] = model[i]
            resid_col[i] = residuals[i]

        # Initialise FITS Table with extension set to obs id
        table = gammalib.GFitsBinTable(nrows)
        # If name not empty add leading blank
        if obs_id != '':
            obs_id = ' ' + obs_id
        table.extname('RESIDUALS ' + obs_id)

        # Add Header card to specify algorithm used for residual computation
        table.card('ALGORITHM', self['algorithm'].string(),
                   'Algorithm used for computation of residuals')

        # Append filled columns to fits table
        table.append(energy_low)
        table.append(energy_high)
        table.append(counts_col)
        table.append(model_col)
        table.append(resid_col)

        return table

    def _append_column(self, table, name, data):
        """
        Append optional column to residual table
        :param table: `~gammalib.GFitsBinTable'
        :param name: str column name
        :param data: `gammalib.GNdarray' table to be filled into new column
        :return: `~gammalib.GFitsBinTable'
        """
        # Check size compatibility
        if table.nrows() == data.size():
            pass
        # Otherwise throw error
        else:
            msg = 'csresspec._append_column: FITS table and data have ' \
                  + 'incompatible size.'
            raise RuntimeError(msg)

        # Create column
        column = gammalib.GFitsTableDoubleCol(name, table.nrows())

        # Fill data
        for i, value in enumerate(data):
            column[i] = value

        # Append new column to table
        table.append(column)

        # Return modified table
        return table

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

        # Stack observations if requested
        if self._stack:
            self.obs(self._stack_observations())

        # Loop over observations and calculate residuals
        for s, observation in enumerate(self.obs()):

            # Retrieve and store obs id
            obs_id = observation.id()
            # If observation id is empty and there is more than one observation
            # replace with incremental number
            if obs_id == '' and self.obs().size() > 1:
                obs_id = str(s)

            # Turn into observation container and assign models
            obs = gammalib.GObservations()
            obs.append(observation)
            obs.models(self.obs().models())

            # if 3D observation
            if obs[0].classname() == 'GCTAObservation':

                ## Prepare Observations

                # If already binned set the evlist_info dictionary to have
                # attribute was_list False
                if obs[0].eventtype() == 'CountsCube':
                    evlist_info = {'was_list': False}
                # Otherwise bin now
                else:
                    # we remember if we binned an event list
                    # so that we can mask only the ROI for residual calculation
                    obs, evlist_info = self._bin_evlist(obs)

                ## Calculate Model and residuals
                # If model cube is provided load it
                if self._use_maps:
                    modcube = gammalib.GCTAEventCube(self['inmodel'].filename())
                # Otherwise calculate it now
                else:
                    modelcube = ctools.ctmodel(obs)
                    modelcube.run()
                    modcube = modelcube.cube().copy()

                # Extract cntcube for residual computation
                cntcube = obs[0].events().copy()

                # Derive count spectra from cubes
                counts = self._cube_to_spectrum(cntcube, evlist_info)
                model = self._cube_to_spectrum(modcube, evlist_info)

                # Calculate residuals
                residuals = obsutils.residuals(self, counts, model)

                # Extract energy bounds
                ebounds = cntcube.ebounds()

                # Fill results table
                table = self._residuals_table(obs_id, ebounds, counts, model,
                                              residuals)

                ## Calculate models of individual components if requested
                if self['components'].boolean():
                    for component in obs.models():
                        # Set model cube models to
                        # individual component
                        model_cont = gammalib.GModels()
                        model_cont.append(component)
                        modelcube.obs().models(model_cont)

                        # Run model cube
                        modelcube.run()

                        # Extract spectrum of individual component
                        modcube = modelcube.cube().copy()
                        model = self._cube_to_spectrum(modcube, evlist_info)

                        # append component to table
                        table = self._append_column(table, component.name(),
                                                    model)

            # otherwise, if On/Off
            elif obs[0].classname() == 'GCTAOnOffObservation':

                ## Calculate Model and residuals

                ## Calculate models of individual components if requested

                pass

            ## Append results table to output file
            self._fits.append(table)

    def save(self):
        """
        Save residuals
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save residuals')

        # Continue only if FITS file is valid
        if self._fits != None:
            # Get outfile parameter
            outfile = self['outfile'].filename()

            # Log file name
            self._log_value(gammalib.NORMAL, 'Residuals file', outfile.url())

            # Save residuals
            self._fits.saveto(outfile, self['clobber'].boolean())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Create instance of application
    app = csresspec(sys.argv)

    # Execute application
    app.execute()
