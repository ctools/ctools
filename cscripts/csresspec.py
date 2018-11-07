#! /usr/bin/env python
# ==========================================================================
# Generates a residual spectrum.
#
# Copyright (C) 2017-2018 Luigi Tibaldo
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
    Generates a residual spectrum
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
        self._stack    = False
        self._mask     = False
        self._fits     = gammalib.GFits()

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
        n_binned   = 0
        n_onoff    = 0
        for obs in self.obs():
            if obs.classname() == 'GCTAObservation':
                if obs.eventtype() == 'CountsCube':
                    n_binned += 1
                else:
                    n_unbinned += 1
            elif obs.classname() == 'GCTAOnOffObservation':
                n_onoff += 1
        n_cta   = n_unbinned + n_binned + n_onoff
        n_other = self.obs().size() - n_cta

        # Query whether to compute model for individual components
        components = self['components'].boolean()

        # If there is only one binned observation and no model for individual
        # components is required, query for precomputed model file and set
        # use_maps to True
        if self.obs().size() == 1 and n_binned == 1 and not components:
            self._use_maps = self['modcube'].is_valid()

        # If there are unbinned observations query the energy binning parameters
        if n_unbinned != 0:
            self['ebinalg'].string()
            if self['ebinalg'].string() == 'FILE':
                self['ebinfile'].filename()
            else:
                self['emin'].real()
                self['emax'].real()
                self['enumbins'].integer()
            if n_cta > n_unbinned:
                n_notunbin = n_cta - n_unbinned

        # If there is more than one observation, and observations are all
        # unbinned or all onoff query user to know if they wish stacked results
        if self.obs().size() > 1 and \
                (n_unbinned == self.obs().size() or n_onoff == self.obs().size()):
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
                self['regfile'].query()

        # Unless all observations are On/Off, or we are using precomputed model
        # maps query whether to use energy dispersion
        if n_onoff == n_cta or self._use_maps:
            pass
        else:
            self['edisp'].boolean()

        # Query algorithm for residual computation
        self['algorithm'].string()

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()

        # Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Write header for observation census
        self._log_header1(gammalib.TERSE, 'Observation census')

        # Log census of input observations
        self._log_value(gammalib.NORMAL, 'Unbinned CTA observations', n_unbinned)
        self._log_value(gammalib.NORMAL, 'Binned CTA observations', n_binned)
        self._log_value(gammalib.NORMAL, 'On/off CTA observations', n_onoff)
        self._log_value(gammalib.NORMAL, 'Other observations', n_other)
        if n_other > 0:
            msg = 'WARNING: Only CTA observation can be handled, all non-CTA ' \
                  + 'observations will be ignored.'
            self._log_string(gammalib.TERSE, msg)

        # Log for unbinned observations
        if n_unbinned != 0:
            msg = ' User defined energy binning will be used for %d unbinned ' \
                  'observations.' % (n_unbinned)
            self._log_string(gammalib.TERSE, msg)
            if n_cta > n_unbinned:
                msg = ' The intrinsic binning will be used for the remaining ' \
                      '%d CTA observations.' % (n_notunbin)
                self._log_string(gammalib.TERSE, msg)

        # Signal how energy dispersion is applied
        if n_onoff == n_cta or self._use_maps:
            msg = ' Energy dispersion is applied based on the input data/model ' \
                  + 'and not according to the edisp parameter'
            self._log_string(gammalib.TERSE, msg)

        # Return
        return

    def _bin_evlist(self, obs):
        """
        Turn single event list into counts cube

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Observation container with single event list

        Returns
        -------
        obs, info : `~gammalib.GObservations`, dict
            Binned observation container and dictionary with event list ROI
            and energy range information
        """
        # Retrieve information about ROI in event list
        roi = obs[0].roi()
        ra  = roi.centre().dir().ra_deg()
        dec = roi.centre().dir().dec_deg()
        rad = roi.radius()

        # We will cover the whole ROI with 0.02 deg binning
        npix = int(2.0 * rad / 0.02) + 1

        # Log binning of events
        self._log_string(gammalib.EXPLICIT, 'Binning events')

        # Bin events
        cntcube = ctools.ctbin(obs)
        cntcube['xref']     = ra
        cntcube['yref']     = dec
        cntcube['binsz']    = 0.02
        cntcube['nxpix']    = npix
        cntcube['nypix']    = npix
        cntcube['proj']     = 'TAN'
        cntcube['coordsys'] = 'CEL'
        cntcube['ebinalg']  = self['ebinalg'].string()
        if self['ebinalg'].string() == 'FILE':
            cntcube['ebinfile'] = self['ebinfile'].filename()
        else:
            cntcube['enumbins'] = self['enumbins'].integer()
            cntcube['emin']     = self['emin'].real()
            cntcube['emax']     = self['emax'].real()
        cntcube.run()

        # Retrieve the binned observation container
        binned_obs = cntcube.obs().copy()

        # Check if energy boundaries provided by user extend beyond the
        # content of the event list
        if self['emin'].real() > obs[0].events().emin().TeV():
            emin = 'INDEF'
        else:
            emin = obs[0].events().emin().TeV()
        if self['emax'].real() < obs[0].events().emax().TeV():
            emax = 'INDEF'
        else:
            emax = obs[0].events().emax().TeV()

        # Log energy range
        self._log_value(gammalib.EXPLICIT, 'Minimum energy (TeV)', emin)
        self._log_value(gammalib.EXPLICIT, 'Maximum energy (TeV)', emax)

        # Put ROI and E bound info in dictionary
        info = {'was_list': True, 'roi_ra': ra, 'roi_dec': dec, 'roi_rad': rad,
                'emin': emin, 'emax': emax}

        # Return new oberservation container
        return binned_obs, info

    def _masked_cube(self, cube, ra, dec, rad, emin='INDEF', emax='INDEF',
                     regfile='NONE'):
        """
        Mask an event cube and returns the masked cube

        Parameters
        ----------
        cube : `~gammalib.GCTAEventCube`
            Event cube
        ra : float (str 'INDEF' for no selection on direction)
            Right Ascension (deg)
        dec : float (str 'INDEF' for no selection on direction)
            Declination (deg)
        rad : float (str 'INDEF' for no selection on direction)
            Radius (deg)
        emin : float (str 'INDEF' for no selection on energy)
            Minimum energy (TeV)
        emax : float (str 'INDEF' for no selection on energy)
            Maximum energy (TeV)

        Returns
        -------
        cube : `~gammalib.GCTAEventCube`
            Event cube
        """
        # Turn cube into observation container to feed to ctcubemask
        obs = gammalib.GCTAObservation()
        obs.events(cube)
        obs_cont = gammalib.GObservations()
        obs_cont.append(obs)

        # Use ctcubemask to mask event cube pixels
        cubemask = ctools.ctcubemask(obs_cont)
        cubemask['ra']      = ra
        cubemask['dec']     = dec
        cubemask['rad']     = rad
        cubemask['emin']    = emin
        cubemask['emax']    = emax
        cubemask['regfile'] = regfile
        cubemask.run()

        # Extract copy of cube from observation container (copy is needed to
        # avoid memory leaks in SWIG)
        cube = cubemask.obs()[0].events().copy()

        # Return cube
        return cube

    def _cube_to_spectrum(self, cube, evlist_info):
        """
        Derive from event cube a count spectrum

        If data come from event list use only the ROI and energy range of
        the original data. Apply user defined mask if requested.

        Parameters
        ----------
        cube : `~gammalib.GCTAEventCube`
            Event cube
        evlist_info : dict
            Dictionary with information on original event list

        Returns
        -------
        array : `~gammalib.GNdarray'
            Counts spectrum
        """
        # If we started from event list mask the ROI only
        # for residual computation
        if evlist_info['was_list']:
            msg = 'Masking ROI from original event list'
            self._log_string(gammalib.EXPLICIT, msg)
            cube = self._masked_cube(cube, evlist_info['roi_ra'],
                                     evlist_info['roi_dec'],
                                     evlist_info['roi_rad'],
                                     emin=evlist_info['emin'],
                                     emax=evlist_info['emax'])

        # Apply user mask
        if self._mask:
            if self['regfile'].is_valid():
                regfile = self['regfile']
            else:
                regfile = 'NONE'
            msg = 'Masking ROI requested by user'
            self._log_string(gammalib.EXPLICIT, msg)
            cube = self._masked_cube(cube, self['ra'], self['dec'],
                                     self['rad'],
                                     regfile=regfile)

        # Extract skymap and clip at 0 to null masked areas
        counts = cube.counts().copy()
        counts = counts.clip(0.)

        # Convert skymap into GNdarray count spectrum
        counts = counts.counts()

        # Return
        return counts

    def _residuals_table(self, obs_id, ebounds, counts, model, residuals):
        """
        Create a Fits Table and store counts, model, and residuals

        Parameters
        ----------
        obs_id : str
            Observation id
        ebounds : `~gammalib.GEbounds'
            Energy boundaries
        counts : `~gammalib.GNdarray'
            Counts spectrum
        model : `~gammalib.GNdarray'
            Model spectrum
        residuals : `~gammalib.GNdarray'
            Residual spectrum

        Returns
        -------
        table : `~gammalib.GFitsBinTable()'
            Residual spectrum as FITS binary table
        """
        # Create FITS table columns
        nrows = ebounds.size()
        energy_low  = gammalib.GFitsTableDoubleCol('Emin', nrows)
        energy_high = gammalib.GFitsTableDoubleCol('Emax', nrows)
        counts_col  = gammalib.GFitsTableDoubleCol('Counts', nrows)
        model_col   = gammalib.GFitsTableDoubleCol('Model', nrows)
        resid_col   = gammalib.GFitsTableDoubleCol('Residuals', nrows)
        energy_low.unit('TeV')
        energy_high.unit('TeV')

        # Fill FITS table columns
        for i in range(nrows):
            energy_low[i]  = ebounds.emin(i).TeV()
            energy_high[i] = ebounds.emax(i).TeV()
            counts_col[i]  = counts[i]
            model_col[i]   = model[i]
            resid_col[i]   = residuals[i]

        # Initialise FITS Table with extension set to obs id
        table = gammalib.GFitsBinTable(nrows)
        table.extname('RESIDUALS' + obs_id)

        # Add Header card to specify algorithm used for residual computation
        table.card('ALGORITHM', self['algorithm'].string(),
                   'Algorithm used for computation of residuals')

        # Append filled columns to fits table
        table.append(energy_low)
        table.append(energy_high)
        table.append(counts_col)
        table.append(model_col)
        table.append(resid_col)

        # Return binary table
        return table

    def _append_column(self, table, name, data):
        """
        Append optional column to residual table

        Parameters
        ----------
        table : `~gammalib.GFitsBinTable'
            FITS binary table
        name : str
            Column name
        data : `~gammalib.GEbounds'
            Data to be filled into new column

        Returns
        -------
        table : `~gammalib.GFitsBinTable'
            FITS binary table
        """
        # Check size compatibility
        if table.nrows() == data.size():
            pass
        # Otherwise throw error
        else:
            msg = 'csresspec._append_column: FITS table and data have ' \
                  'incompatible size.'
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

    def _residuals_3D(self, obs, obs_id, ccube='NONE'):
        """
        Calculate residuals for 3D observation

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Observation container with a single observations of type GCTAObservation
        obs_id : str
            Observation ID
        ccube : `~gammalib.GCTAEventCube'
            Count cube with stacked events lists

        Returns
        -------
        table : `~gammalib.GFitsBinTable()'
            Residual spectrum as FITS binary table
        """

        # If binned data already exist set the evlist_info dictionary to have
        # attribute was_list False
        if obs[0].eventtype() == 'CountsCube' or ccube!='NONE':
            evlist_info = {'was_list': False}

        # ... otherwise bin now
        else:
            # we remember if we binned an event list so that we can
            # mask only the ROI for residual calculation
            msg = 'Setting up binned observation'
            self._log_string(gammalib.NORMAL, msg)
            obs, evlist_info = self._bin_evlist(obs)

        # Calculate Model and residuals. If model cube is provided load
        # it
        if self._use_maps:
            modcube = gammalib.GCTAEventCube(self['modcube'].filename())

        # ... otherwise calculate it now
        else:
            msg = 'Computing model cube'
            self._log_string(gammalib.NORMAL, msg)
            modelcube = ctools.ctmodel(obs)
            if ccube != 'NONE':
                modelcube.cube(ccube)
            modelcube['edisp'] = self['edisp'].boolean()
            modelcube.run()
            modcube = modelcube.cube().copy()

        # Extract cntcube for residual computation
        if ccube != 'NONE':
            cntcube = ccube
        else:
            cntcube = obs[0].events().copy()

        # Derive count spectra from cubes
        msg = 'Computing counts, model, and residual spectra'
        self._log_string(gammalib.NORMAL, msg)
        counts = self._cube_to_spectrum(cntcube, evlist_info)
        model  = self._cube_to_spectrum(modcube, evlist_info)

        # Calculate residuals
        residuals = obsutils.residuals(self, counts, model)

        # Extract energy bounds
        ebounds = cntcube.ebounds()

        # Fill results table
        msg = 'Filling residual table'
        self._log_string(gammalib.NORMAL, msg)
        table = self._residuals_table(obs_id, ebounds, counts, model,
                                      residuals)

        # Calculate models of individual components if requested
        if self['components'].boolean():
            for component in obs.models():
                # Log component
                self._log_value(gammalib.NORMAL,
                                'Computing model component',
                                component.name())

                # Set model cube models to individual component
                model_cont = gammalib.GModels()
                model_cont.append(component)
                modelcube.obs().models(model_cont)

                # Reset base cube that was modified internally by ctmodel
                if ccube != 'NONE':
                    modelcube.cube(ccube)

                # Run model cube
                modelcube['edisp'] = self['edisp'].boolean()
                modelcube.run()

                # Extract spectrum of individual component
                modcube = modelcube.cube().copy()
                model = self._cube_to_spectrum(modcube, evlist_info)

                # append component to table
                table = self._append_column(table, component.name(),
                                            model)

        return table

    def _residuals_OnOff(self,obs,obs_id):
        """
        Calculate residual for OnOff observation

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Observation container with a single observations of type GCTAOnOffObservation
        obs_id : str
            Observation ID

        Returns
        -------
        table : `~gammalib.GFitsBinTable'
            Residual spectrum as FITS binary table
        """

        # Calculate Counts, Model and residuals
        msg = 'Computing counts, model, and residual spectra'
        self._log_string(gammalib.NORMAL, msg)

        onoff = obs[0]

        # On spectrum
        counts = onoff.on_spec().counts_spectrum()

        # Model for On region
        background = onoff.model_background(obs.models()).counts_spectrum()
        alpha      = onoff.on_spec().backscal_spectrum()
        model      = background.copy()
        model     *= alpha
        model     += onoff.model_gamma(obs.models()).counts_spectrum()

        # On Residuals
        residuals = obsutils.residuals(self, counts, model)

        # Extract energy bounds
        ebounds = onoff.on_spec().ebounds()

        # Fill results table
        msg = 'Filling residual table'
        self._log_string(gammalib.NORMAL, msg)
        table = self._residuals_table(obs_id, ebounds, counts, model,
                                      residuals)

        # Get Off spectrum and add to table
        msg = 'Computing counts, model, and residual spectra for Off regions'
        self._log_string(gammalib.NORMAL, msg)
        counts_off = onoff.off_spec().counts_spectrum()
        table = self._append_column(table, 'Counts_Off',
                                    counts_off)

        # Add background/Off model to table
        table = self._append_column(table, 'Model_Off',
                                    background)

        # Calculate Off residuals and add to table
        residuals_off = obsutils.residuals(self, counts_off, background)
        table = self._append_column(table, 'Residuals_Off',
                                    residuals_off)

        # Calculate models of individual components if requested
        if self['components'].boolean():
            for component in obs.models():
                # If background pass
                # We always add the background at the end so that
                # we accommodate WSTAT for which the background is not
                # mandatory in the model
                if component.classname() == 'GCTAModelIrfBackground':
                    pass
                # Otherwise calculate gamma component and append to Table
                else:
                    self._log_value(gammalib.NORMAL,
                                    'Computing model for component',
                                    component.name())
                    # Create observation container for individual components
                    model_cont = gammalib.GModels()
                    model_cont.append(component)
                    # Calculate gamma model
                    model = onoff.model_gamma(model_cont)
                    model = model.counts_spectrum()
                    # Append to table
                    table = self._append_column(table, component.name(),
                                                model)
            # Add background already calculated
            self._log_value(gammalib.NORMAL,
                            'Computing model for component',
                            'Background')
            background *= onoff.on_spec().backscal_spectrum()
            table = self._append_column(table, 'Background',
                                        background)

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

        # Stack On/Off observations if requested
        if self._stack and self.obs()[0].classname() == 'GCTAOnOffObservation':
            msg = 'Stacking %d On/Off observations.' % (self.obs().size())
            self._log_string(gammalib.NORMAL, msg)
            stacked = gammalib.GCTAOnOffObservation(self.obs())
            stacked_obs = gammalib.GObservations()
            stacked_obs.append(stacked)
            stacked_obs.models(self.obs().models())
            self.obs(stacked_obs)

        # Log processing header
        self._log_header1(gammalib.TERSE, 'Processing Observations')

        #If observations are unbinned and we stack
        if self._stack and self.obs()[0].classname() == 'GCTAObservation':

            msg = 'Computing count cube from multiple unbinned observations'
            self._log_string(gammalib.NORMAL, msg)

            # Build count cube
            binning = ctools.ctbin(self.obs())
            binning['xref']     = self['xref'].real()
            binning['yref']     = self['yref'].real()
            binning['proj']     = self['proj'].string()
            binning['coordsys'] = self['coordsys'].string()
            binning['ebinalg']  = self['ebinalg'].string()
            binning['nxpix']    = self['nxpix'].integer()
            binning['nypix']    = self['nypix'].integer()
            binning['binsz']    = self['binsz'].real()
            if self['ebinalg'].string() == 'FILE':
                binning['ebinfile'] = self['ebinfile'].filename().file()
            else:
                binning['enumbins'] = self['enumbins'].integer()
                binning['emin']     = self['emin'].real()
                binning['emax']     = self['emax'].real()
            binning['chatter'] = self['chatter'].integer()
            binning['clobber'] = self['clobber'].boolean()
            binning['debug']   = self['debug'].boolean()
            binning.run()

            #compute residuals using cube
            table = self._residuals_3D(self.obs(),'',binning.cube())

            # Append results table to output file
            self._fits.append(table)

        # Otherwise, loop over observations and calculate residuals
        else:
            for s, observation in enumerate(self.obs()):

                # Retrieve and store obs id
                obs_id = observation.id()

                # If observation id is empty and there is more than one observation
                # replace with incremental number
                if obs_id == '' and self.obs().size() > 1:
                    obs_id = str(s)

                # Log processing of observation
                if self.obs().size() > 1:
                    self._log_header2(gammalib.NORMAL,
                                      'Processing observation %s' % obs_id)

                # Turn into observation container and assign models
                obs = gammalib.GObservations()
                obs.append(observation)
                obs.models(self.obs().models())

                # If 3D observation
                if obs[0].classname() == 'GCTAObservation':
                    table = self._residuals_3D(obs,obs_id)

                # otherwise, if On/Off
                elif obs[0].classname() == 'GCTAOnOffObservation':
                    table = self._residuals_OnOff(obs,obs_id)

                # Append results table to output file
                self._fits.append(table)

        # Return
        return

    def save(self):
        """
        Save residuals
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save residuals')

        # Get outfile parameter
        outfile = self['outfile'].filename()

        # Log file name
        self._log_value(gammalib.NORMAL, 'Residuals file', outfile.url())

        # Save residuals
        self._fits.saveto(outfile, self['clobber'].boolean())

        # Return
        return

    def resspec(self):
        """
        Return residual FITS file

        Returns
        -------
        fits : `~gammalib.GFits'
            FITS file containing residuals
        """
        # Return
        return self._fits


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csresspec(sys.argv)

    # Execute application
    app.execute()
