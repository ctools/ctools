#!/usr/bin/env python
# ==========================================================================
# Selects observations from an observation definition XML file
#
# Copyright (C) 2016 Juergen Knoedlseder
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


# ================= #
# csobsselect class #
# ================= #
class csobsselect(ctools.cscript):
    """
    Selects observations from an observation definition XML file
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'csobsselect'
        self._version = '1.2.0'

        # Initialise observation container from constructor arguments
        self._obs, argv = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup observations
        """
        # Query input parameters
        self['inobs'].filename()
        
        # Query relevant pointing selection parameters
        pntselect = self['pntselect'].string()
        coordsys  = self['coordsys'].string()
        if coordsys == 'CEL':
            self['ra'].real()
            self['dec'].real()
        else:
            self['glon'].real()
            self['glat'].real()
        if pntselect == 'CIRCLE':
            self['rad'].real()
        else:
            self['width'].real()
            self['height'].real()

        # Query ahead output model filename
        if self._read_ahead():
            self['outobs'].filename()

        # If there are no observations in container then get them from the
        # parameter file
        if self._obs.size() == 0:
            self._obs = self._get_observations(False)

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _select_observation(self, obs):
        """
        Select observation

        Parameters
        ----------
        obs : `~gammalib.GObservation`
            Observation

        Returns
        -------
        select : bool
            Observation selection flag
        """
        # Initialise selection flag
        select = False

        # Continue only if observation is a CTA observation
        if obs.classname() == 'GCTAObservation':

            # Get selection type
            pntselect = self['pntselect'].string()

            # Dispatch according to selection type
            if pntselect == 'CIRCLE':
                select = self._select_circle(obs)
            else:
                select = self._select_box(obs)

        # Return selection flag
        return select

    def _select_circle(self, obs):
        """
        Select observation in a pointing circle

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation`
            CTA observation

        Returns
        -------
        select : bool
            Observation selection flag
        """
        # Get coordinate system of selection circle
        coordsys = self['coordsys'].string()

        # Get pointing direction
        pnt = obs.pointing().dir()

        # Set selection circle centre
        centre = gammalib.GSkyDir()
        if coordsys == 'CEL':
            centre.radec_deg(self['ra'].real(), self['dec'].real())
        else:
            centre.lb_deg(self['glon'].real(), self['glat'].real())

        # Set selection flag according to distance
        if centre.dist_deg(pnt) <= self['rad'].real():
            select = True
            msg    = 'inside selection circle'
        else:
            select = False
            msg    = 'outside selection circle'

        # Set observation name
        if len(obs.name()) > 0:
            name = obs.name()
        else:
            name = 'Observation'

        # Prepend obervation ID
        if len(obs.id()) > 0:
            msg = obs.id()+' '+msg

        # Log observation selection
        self._log_value(gammalib.NORMAL, name, msg)

        # Return selection flag
        return select

    def _select_box(self, obs):
        """
        Select observation in a pointing box

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation`
            CTA observation

        Returns
        -------
        select : bool
            Observation selection flag
        """
        # Initialise selection flag
        select = False
        msg    = 'outside selection box'

        # Get coordinate system of selection circle
        coordsys = self['coordsys'].string()

        # Get pointing direction
        pnt = obs.pointing().dir()

        # Get selection box width and height
        width  = self['width'].real()
        height = self['height'].real()

        # Make selection for celestial coordinates ...
        if coordsys == 'CEL':

            # Determine box boundaries
            ra      = self['ra'].real()
            dec     = self['dec'].real()
            ra_min  = ra - 0.5 * width
            ra_max  = ra + 0.5 * width
            dec_min = dec - 0.5 * height
            dec_max = dec + 0.5 * height

            # Check if pointing lies within boundaries
            ra       = pnt.ra_deg()
            ra_plus  = ra + 360.0
            ra_minus = ra - 360.0
            dec = pnt.dec_deg()
            if dec >= dec_min and dec <= dec_max:
                if ra       >= ra_min and ra       <= ra_max or \
                   ra_plus  >= ra_min and ra_plus  <= ra_max or \
                   ra_minus >= ra_min and ra_minus <= ra_max:
                    select = True
                    msg    = 'inside selection box'

        # ... or galactic coordinates
        else:

            # Determine box boundaries
            glon     = self['glon'].real()
            glat     = self['glat'].real()
            glon_min = glon - 0.5 * width
            glon_max = glon + 0.5 * width
            glat_min = glat - 0.5 * height
            glat_max = glat + 0.5 * height

            # Check if pointing lies within boundaries
            glon       = pnt.l_deg()
            glon_plus  = glon + 360.0
            glon_minus = glon - 360.0
            glat       = pnt.b_deg()
            if glat >= glat_min and glat <= glat_max:
                if glon       >= glon_min and glon       <= glon_max or \
                   glon_plus  >= glon_min and glon_plus  <= glon_max or \
                   glon_minus >= glon_min and glon_minus <= glon_max:
                    select = True
                    msg    = 'inside selection box'

        # Set observation name
        if len(obs.name()) > 0:
            name = obs.name()
        else:
            name = 'Observation'

        # Prepend obervation ID
        if len(obs.id()) > 0:
            msg = obs.id()+' '+msg

        # Log observation selection
        self._log_value(gammalib.NORMAL, name, msg)

        # Return selection flag
        return select


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

        # Initialise empty observation container for selected observations
        selected_obs = gammalib.GObservations()

        # Write input observation container into logger
        self._log_observations(gammalib.NORMAL, self._obs, 'Input observation')

        # Write header
        self._log_header1(gammalib.NORMAL, 'Observation selection')

        # Loop over observations
        for obs in self._obs:

            # If observation is selected then append observation
            if self._select_observation(obs):
                selected_obs.append(obs)

        # Copy selected observations into observation
        self._obs = selected_obs

        # Write input observation container into logger
        self._log_observations(gammalib.NORMAL, self._obs, 'Selected observation')

        # Return
        return

    def save(self):
        """ 
        Save observation definition XML file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save observations')

        # Get output filename
        outobs = self['outobs'].filename()

        # If file exists and clobber flag is false then raise an exception
        if outobs.exists() and not self._clobber:
            msg = ('Cannot save "'+outobs.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:
            # Log filename
            self._log_value(gammalib.NORMAL, 'Observation definition XML file',
                                             outobs.url())

            # Save observations
            self._obs.save(outobs)

        # Return
        return

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save observation definition file
        self.save()

        # Return
        return    


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csobsselect(sys.argv)

    # Execute application
    app.execute()
