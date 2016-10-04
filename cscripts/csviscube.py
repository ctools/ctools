#! /usr/bin/env python
# ==========================================================================
# Compute a visibility cube
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
import glob
import os
import math
import gammalib
import ctools


# =============== #
# csviscube class #
# =============== #
class csviscube(ctools.cscript):
    """
    Compute a visibility cube
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        
        Parameters
        ----------
        argv : list of str
            List of IRAF command line parameter strings of the form
            ``parameter=3``.

        Raises
        ------
        TypeError
            An invalid number of command line arguments was provided.
        """
        # Set name
        self._name    = 'csviscube'
        self._version = '1.0.0'

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get all parameters
        """
        # Query parameters
        self['tmin'].real()
        self['tmax'].real()
        self['geolon'].real()
        self['geolat'].real()
        self['sunzenith'].real()
        self['moonzenith'].real()
        self['outfile'].filename()
        self['binsz'].real()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _hour_angle_weight(self):
        """
        Compute hour angle weight vector
        
        Compute during which time the array was observing due to hour angle
        constraints. We take here a very simple model where the night lasts
        between 6 and 10 hours (or 90-150 degrees), with a sinusoidal variation
        along the year. Every day the hour angle of the Sun set is advancing by
        4 minutes.
        """
        # Set number of hour angle bins
        nsteps = 1000
        d_min  =  6.0/24.0*nsteps
        d_max  = 10.0/24.0*nsteps
    
        # Initialise hour angle array and weight
        hours  = [0.0 for i in range(nsteps)]
        weight = 24.0/float(nsteps)  # Weight per hour angle step

        # Loop over all days in one year
        for day in range(nsteps):
            h_set  = day
            h_rise = day + \
                     int(d_min + (d_max-d_min)*math.sin(day/float(nsteps)*gammalib.pi))
            for i in range(h_set,h_rise):
                if i >= nsteps:
                    hours[i-nsteps] += weight
                else:
                    hours[i] += weight

        # Return
        return hours

    def _visibility_cube(self):
        """
        Compute visibility cube
        
        Compute visibility cube by displacing the zenith angle map for all
        hour angles. The visibility cube contains the number of hours a given
        celestial position is visible under a given zenith angle interval.
        Summing over all zenith angle intervals specifies for how long a given
        celestial position will be visible.
        """
        # Compute zenith angle map
        map = self._zenith_angle_map()
        
        # Set zenith angle range and bin size
        dz = 1.0
        nz = int(60.0/dz)
    
        # Initialise visibility cube
        self._cube = gammalib.GSkyMap('CAR','CEL',0.0,0.0,-1.0,1.0,360,180,nz)

        # Setup hour angle weights. They specify for how many hours a given
        # hour angle will be observed during one year.
        hours = self._hour_angle_weight()

        # Compute normalisation factor
        norm = 365.0/float(len(hours))

        # Loop over all hour angle weights
        for h, weight in enumerate(hours):

            # Compute temporal shift
            shift = 360.0/float(len(hours)) * float(h)

            # Initialise shifted zenith angle map
            map_shift = gammalib.GSkyMap('CAR','CEL',shift,0.0,-1.0,1.0,360,180)

            # Merge zenith angle map into shifted map
            map_shift += map

            # Loop over all pixels of shifted map
            for i, zenith in enumerate(map_shift):
                iz     = int(zenith/dz)
                if iz < nz:
                    self._cube[i,iz] += weight * norm

        # Return
        return

    def _zenith_angle_map(self):
        """
        Compute zenith angle map

        The zenith angle of a position (ra,dec) depends on the declination and
        the hour angle h and is given by

        zenith(h,dec) = arccos( sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(h) )

        The hour angle h (or local hour angle, LHA) is defined as the difference
        between local siderial time (LST) and the Right Ascension
    
        h = LST - ra

        LST is the difference between the Greenwich siderial time (GST) and the
        geographic longitude of the observer
    
        LST = GST - lon
    
        where lon is counted positive west from the meridian. Hence
    
        h = GST - lon - ra

        The map is computed for h=-ra which is equivalent to GST=lon (or LST=0).
        In other words, the map corresponds to the time when ra=0 passes the
        local meridian.
        """
        # Initialise zenith angle map
        map = gammalib.GSkyMap('CAR','CEL',0.0,0.0,-1.0,1.0,360,180)

        # Set hour angle and declination vectors
        hours = [float(i) for i in range(360)]
        decs  = [float(i)-89.5 for i in range(180)]

        # Get array geographic longitude and latitude
        geolon = self['geolon'].real()
        geolat = self['geolat'].real()

        # Precompute latitude terms
        cos_lat = math.cos(geolat*gammalib.deg2rad)
        sin_lat = math.sin(geolat*gammalib.deg2rad)

        # Loop over all declination and hour angles and compute the zenith
        # angle. Store the zenith angle in the map
        index = 0
        for dec in decs:
            cos_dec = math.cos(dec*gammalib.deg2rad)
            sin_dec = math.sin(dec*gammalib.deg2rad)
            for h in hours:
                cos_h  = math.cos(h*gammalib.deg2rad)
                zenith = math.acos(sin_lat*sin_dec +
                                   cos_lat*cos_dec*cos_h)*gammalib.rad2deg
                map[index] = zenith
                index     += 1

        # Return zenith angle map
        return map


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

        # Compute visibility cube
        self._visibility_cube()

        # Return
        return

    def save(self):
        """
        Save the visibility cube
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save visibility cube')

        # Get outfile parameter
        outfile = self['outfile'].filename()
        
        # Log file name
        self._log_value(gammalib.TERSE, 'Visibility cube file', outfile.url())

        # Save the visibility cube
        self._cube.save(outfile, self['clobber'].boolean())

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

        # Save the cube
        self.save()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csviscube(sys.argv)

    # Execute application
    app.execute()
