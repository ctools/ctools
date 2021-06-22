#! /usr/bin/env python
# ==========================================================================
# Select observations from COMPTEL database
#
# Copyright (C) 2019-2021 Juergen Knoedlseder
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
import glob
import os
import gammalib
import ctools


# ================== #
# comobsselect class #
# ================== #
class comobsselect(ctools.csobservation):
    """
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
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
        if self['tmin'].is_valid() and self['tmax'].is_valid():
            self['tmin'].time()
            self['tmax'].time()

        # Query ahead output model filename
        if self._read_ahead():
            self['outobs'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _load_dbase(self):
        """
        Load database
        """
        # Get database FITS filename
        filename = self['dbase'].filename()
    
        # Load database FITS file
        fits = gammalib.GFits(filename)
        
        # Get FITS table
        table = fits.table('XML')
    
        # Return FITS table (need a copy here because FITS file goes out of
        # scope after leaving the method)
        return table.copy()

    def _select_observation(self, pnt, tstart, tstop):
        """
        Select observation

        Parameters
        ----------
        pnt : `~gammalib.GSkyDir`
            Observation pointing direction
        tstart : `~gammalib.GTime`
            Observation start time
        tstop : `~gammalib.GTime`
            Observation stop time

        Returns
        -------
        select, msg : bool, string
            Observation selection flag and message
        """
        # Initialise selection flag and message
        select = True
        msg    = ''

        # If there is a valid User start and stop time and the observation
        # stop is before the start time or the observation start after the
        # stop time, then skip the observation
        if self['tmin'].is_valid() and self['tmax'].is_valid():
            tmin = self['tmin'].time()
            tmax = self['tmax'].time()
            if tstop < tmin or tstart > tmax:
                msg    = 'outside time interval'
                select = False

        # If temporally selection was passed then now select spatially
        if select:

            # Get selection type
            pntselect = self['pntselect'].string()

            # Dispatch according to selection type
            if pntselect == 'CIRCLE':
                select, msg = self._select_circle(pnt)
            else:
                select, msg = self._select_box(pnt)

        # Return selection flag and message
        return select, msg

    def _select_circle(self, pnt):
        """
        Select observation in a pointing circle

        Parameters
        ----------
        pnt : `~gammalib.GSkyDir`
            Pointing direction

        Returns
        -------
        select, msg : bool, string
            Observation selection flag and message
        """
        # Get coordinate system of selection circle
        coordsys = self['coordsys'].string()

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

        # Return selection flag and message
        return select, msg

    def _select_box(self, pnt):
        """
        Select observation in a pointing box

        Parameters
        ----------
        pnt : `~gammalib.GSkyDir`
            Pointing direction

        Returns
        -------
        select, mst : bool, string
            Observation selection flag and message
        """
        # Initialise selection flag
        select = False
        msg    = 'outside selection box'

        # Get coordinate system of selection circle
        coordsys = self['coordsys'].string()

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

        # Return selection flag and message
        return select, msg


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

        # Clear selected observations
        self.obs().clear()

        # Load database
        dbase = self._load_dbase()

        # Loop over all entries
        for i in range(dbase.nrows()):

            # Set VP name
            name = dbase['VP'].string(i)

            # Initialise message string
            msg = ''

            # Get pointing direction
            glon = dbase['GLON_SCZ'].real(i)
            glat = dbase['GLAT_SCZ'].real(i)
            pnt  = gammalib.GSkyDir()
            pnt.lb_deg(glon, glat)

            # Get start and stop time
            tjd_start  = dbase['TJD_START'].integer(i)
            tics_start = dbase['TICS_START'].integer(i)
            tjd_stop   = dbase['TJD_STOP'].integer(i)
            tics_stop  = dbase['TICS_STOP'].integer(i)
            tstart     = gammalib.com_time(tjd_start, tics_start)
            tstop      = gammalib.com_time(tjd_stop,  tics_stop)

            # Check for selection
            select, msg = self._select_observation(pnt, tstart, tstop)

            # Append observation if selected
            if select:

                # Load observation
                filename = dbase['XML'].string(i)
                obs      = gammalib.GObservations(filename)

                # Append observation
                self.obs().append(obs[0])

            # Log observation selection
            self._log_value(gammalib.NORMAL, name, msg)

        # Log selected observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

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
        if outobs.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outobs.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:
            # Log filename
            self._log_value(gammalib.NORMAL, 'Obs. definition XML file',
                                             outobs.url())

            # Save observations
            self.obs().save(outobs)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comobsselect(sys.argv)

    # Execute application
    app.execute()
