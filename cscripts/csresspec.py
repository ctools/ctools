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
        self._components = False
        self._use_maps = False
        self._stack = False

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Set observation if not done before
        if self.obs().is_empty():
            self._require_inobs('csresspec._get_parameters()')
            self.obs(self._get_observations())

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

        # Query wether to compute model for individual components
        self._components = self['components'].boolean()

        # If there is only one binned observation and no model for individual
        # components is required, query for precomputed model file and set
        # use_maps to True
        if self.obs().size() == 1 and n_binned == 1 and not self._components:
            modcube = self['modcube'].filename()
            if modcube != 'NONE':
                self._use_maps = True

        # If we are not using a precomputed model query input XML model file
        if not self._use_maps:
            modxml = self['inmodel'].filename()
            # if None check whether models are provided in observation
            # container, otherwise throw exception and stop here
            if modxml == 'NONE' and self.obs().models().is_empty():
                msg = 'No model provided. Please specify an input XML model file'
                raise RuntimeError(msg)

        # If there are unbinned observations query the energy binning parameters
        if n_unbinned != 0:
            self['ebinalg'].string()
            if self['ebinalg'].string() == 'FILE':
                self['ebinfile'].filename()
            else:
                self['emin'].real()
                self['emax'].real()
                self['enumbins'].integer()

        # If there is more than one observation, and observations are unbinned
        # or onoff query user to know if they wish stacked results
        if self.obs().size() > 1 and (
                        n_unbinned == self.obs().size() or n_onoff == self.obs().size()):
            self['stack'].boolean()

        # Apply energy dispersion
        self['edisp'].boolean()

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()

        # Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Log census of input observations
        self._log_value(gammalib.TERSE, 'Unbinned CTA observations', n_unbinned)
        self._log_value(gammalib.TERSE, 'Binned CTA observations', n_binned)
        self._log_value(gammalib.TERSE, 'On/off CTA observations', n_onoff)
        self._log_value(gammalib.TERSE, 'Other observations', n_other)

        # Return
        return
