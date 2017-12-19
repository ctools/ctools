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
import math


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
            self._log_value(gammalib.TERSE, msg)
            if n_cta > n_unbinned:
                n_notunbin = n_cta - n_unbinned
                msg = 'The intrinsic binning will be used for the remaining %s CTA observations.' % (
                    n_notunbin)
                self._log_value(gammalib.TERSE, msg)

        # If there is more than one observation, and observations are all
        # unbinned or all onoff query user to know if they wish stacked results
        if self.obs().size() > 1 and (
                        n_unbinned == self.obs().size() or n_onoff == self.obs().size()):
            self._stack = self['stack'].boolean()
            # If we are to stack event lists query parameters for cube creation
            if n_unbinned == self.obs().size():
                self['coordsys'].string()
                self['proj'].string()
                self['xref'].real()
                self['yref'].real()
                self['nxpix'].real()
                self['nypix'].real()
                self['binsz'].real()

        # If we are not using a precomputed model query input XML model file
        if not self._use_maps:
            modxml = self['inmodel'].filename()
            # If None check whether models are provided in observation
            # container, otherwise throw exception and stop here
            if modxml == 'NONE' and self.obs().models().is_empty():
                msg = 'No model provided. Please specify an input XML model file'
                raise RuntimeError(msg)
            # Otherwise set the input model
            else:
                self.obs().models(modxml)

        # Unless all observations are On/Off query for mask definition
        if n_onoff == n_cta:
            pass
        else:
            self._mask = self['mask'].boolean()
            if self._mask:
                self['ra'].real()
                self['dec'].real()
                self['rad'].real()

        # Unless all observations are On/Off, or we are using precomputed model
        # maps query whether to use energy dispersion
        if n_onoff == n_cta or self._use_maps:
            msg = 'Energy dispersion is applied based on the input data/model ' \
                  + 'and not according to the edisp parameter'
            self._log_value(gammalib.TERSE, msg)
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

        :return: stacked_obs container with a single stacked observation
        """
        # Use first observation to determine the type and apply the
        # appropriate stacking method

        # If On/Off use the GCTAOnOffObservation constructor
        if self.obs()[0].classname() == 'GCTAOnOffObservation':
            msg = 'Stacking %s On/Off observations.' % (self.obs().size())
            self._log_value(gammalib.TERSE, msg)
            stacked_obs = gammalib.GCTAOnOffObservation(self.obs())

        # If event list then bin observations
        elif self.obs()[0].classname() == 'GCTAObservation':
            msg = 'Stacking %s event lists.' % (self.obs().size())
            self._log_value(gammalib.TERSE, msg)
            stacked_obs = obsutils.get_stacked_obs(self, self.obs())

        # Return
        return stacked_obs

    def _bin_evlist(self, obs):
        """
        Turn single event list into counts cube
        :param   obs: observation container with single event list
        :return: binned_obs: binned observation container
                 ra: RA of ROI centre
                 dec: Dec of ROI centre
                 rad: ROI radius
                 emin: minimum energy of events
                 emax: maximum energy of events

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
        if self['ebinalg'].string() == 'FILE'
            cntcube['ebinfile'] = self['ebinfile'].filename()
        else:
            cntcube['enumbins'] = self['enumbins'].integer()
            cntcube['emin'] = self['emin'].real()
            cntcube['emax'] = self['emax'].real()
        cntcube.run()

        # Retrieve a new oberservation container
        binned_obs = cntcube.obs().copy()

        # Retrieve models
        models = obs.models()

        # Set models for new oberservation container
        binned_obs.models(models)

        # Check if energy boundaries provided by user extend beyond
        # the content of the event list
        if self['emin'].real() > obs[0].emin().TeV:
            emin = 'INDEF'
        else:
            emin = obs[0].emin().TeV
        if self['emax'].real() < obs[0].emax().TeV:
            emax = 'INDEF'
        else:
            emax = obs[0].emax().TeV

        # Return new oberservation container
        return binned_obs, ra, dec, rad, emin, emax

    def _masked_cube(self, cube, ra, dec, rad, emin='INDEF', emax='INDEF',
                     regfile='NONE'):
        """
        Mask an event cube and returns the masked cube
        :param cube: Events cube
        :param ra:
        :param dec:
        :param rad:
        :param emin:
        :param emax:
        :return: outcube: masked cube
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

        return cubemask.cube()


def _add_binned_response(self, obs):
    """
    Add response to binned observation
    :param obs: observation container
    :return: obs: observation container with response
    """

    # use counts cube to extract geometry
    cube = obs[0].events()
    proj = cube.counts().projection()

    response = obsutils.get_stacked_response(obs,
                                             proj.crval(0),
                                             proj.crval(1),
                                             binsz=max(abs(proj.cdelt(0)),
                                                       abs(proj.cdelt(1))),
                                             # always ensures coverage of whole cube
                                             nxpix=cube.nx(),
                                             nypix=cube.ny(),
                                             coordsys=proj.coordsys(),
                                             proj=proj.code(),
                                             emin=cube.emin().TeV(),
                                             emax=cube.emax().TeV(),
                                             edisp=self['edisp'].boolean())

    if self['edisp'].boolean():
        obs[0].response(response['expcube'], response['psfcube'],
                        response['edispcube'], response['bkgcube'])
    else:
        obs[0].response(response['expcube'], response['psfcube'],
                        response['bkgcube'])

    # Get new models
    models = response['models']

    # Set models for observation container
    obs.models(models)

    # Return observation container
    return obs

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
    for obs in self.obs():

        # Turn into observation container and assign models
        obs = gammalib.GObservations(obs)
        obs.models(self.obs().models())

        # if 3D observations
        if obs[0].classname() == 'GCTAObservation':

            ## Prepare Observations

            # If already binned pass
            was_list = False
            if obs[0].eventtype() == 'CountsCube':
                pass
            # Otherwise bin now
            else:
                # we remember if we binned an event list
                # so that we can mask only the ROI for residual calculation
                was_list = True
                obs, ev_ra, ev_dec, ev_rad, ev_emin, ev_emax = self._bin_evlist(
                    obs)

            # If binned response is present or model cube provided pass
            if obs[0].has_response() or self._use_maps:
                pass
            # Otherwise calculate binned response now
            else:
                obs = self._add_binned_response(obs)

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

            # If we started from event list mask the ROI only
            # for residual computation
            if was_list:
                cntcube = self._masked_cube(cntcube, ev_ra, ev_dec, ev_rad,
                                            emin=ev_emin, emax=ev_emax)
                modcube = self._masked_cube(modcube, ev_ra, ev_dec, ev_rad,
                                            emin=ev_emin, emax=ev_emax)
            else:
                pass

            # Apply user mask
            if self._mask:
                cntcube = self._masked_cube(cntcube, self['ra'], self['dec'],
                                            self['rad'],
                                            regfile=self['regfile'])
                modcube = self._masked_cube(modcube, self['ra'], self['dec'],
                                            self['rad'],
                                            regfile=self['regfile'])
            else:
                pass

            # Calculate residuals

            # Arrays with counts and model
            counts = cntcube.counts().total_counts()
            model =  modcube.counts().total_counts()

            # Get residual algorithm type
            algorithm = self['algorithm'].string()

            # Subtract
            if algorithm == 'SUB':
                residuals = counts - model

            # Subtract and divide by model map
            elif algorithm == 'SUBDIV':
                residuals = counts - model
                residuals /= model

            # Subtract and divide by sqrt of model map
            elif algorithm == 'SUBDIVSQRT':
                residuals = counts - model
                residuals /= model.sqrt()

            # Calculate significance from Li&Ma derivation
            elif algorithm == 'SIGNIFICANCE':
                residuals = counts.copy()

                # Compute sign array
                sign = (counts - model).sign() ##not implemented!!

                # Loop over energy bins
                for i in range(counts.size()):

                    # If the model value > 0.0 do the computation as normal ...
                    model_val = model[i]
                    if model_val > 0.0:

                        # If the data value is also > 0 then compute the
                        # significance^2 and fill it in the residuals ...
                        data_val = counts[i]
                        if data_val > 0.0:
                            log_val = math.log(data_val / model_val)
                            residuals[i] = (
                                              data_val * log_val) + model_val - data_val

                        # ... otherwise compute the reduced value of the above
                        # expression. This is necessary to avoid computing log(0).
                        else:
                            residuals[i] = model_val

                    # ... otherwise hard-code the significance to 0
                    else:
                        residuals[i] = 0.0

                # Compute significance
                residuals *= 2.0
                residuals = residuals.sqrt()
                residuals = residuals * sign

            ## Calculate models of individual components if requested

        # otherwise, if On/Off
        elif obs[0].classname() == 'GCTAOnOffObservation':

            ## Calculate Model and residuals

            ## Calculate models of individual components if requested

            pass
