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
        Constructor
        
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

        # Initialise members
        self._cube = gammalib.GSkyMap()

        # Initialise application by calling the appropriate class constructor
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
        self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _zenith_angle_map(self, nx, ny, dx, dy):
        """
        Compute zenith angle map

        Parameters
        ----------
        nx : int
            Number of Right Ascension pixels.
        ny : int
            Number of Declination pixels.
        dx : float
            Right Ascension pixel size in degrees.
        dy : float
            Declination pixel size in degrees.

        Returns
        -------
        zmap : `~gammalib.GSkyMap`
            Allsky map comprising the zenith angle for an hour angle of 0.

        The zenith angle of a position (ra,dec) depends on the declination and
        the hour angle h and is given by

        zenith(h,dec) = arccos( sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(h) )

        The hour angle h (or local hour angle, LHA) is defined as the difference
        between local siderial time (LST) and the Right Ascension
    
        h = LST - ra

        The map is computed for h=-ra which is equivalent to GST=lon (or LST=0).
        In other words, the map corresponds to the time when ra=0 passes through
        the local meridian.
        """
        # Initialise zenith angle map
        zmap = gammalib.GSkyMap('CAR','CEL',0.0,0.0,-dx,dy,nx,ny)

        # Set hour angle and declination vectors
        hours = [(float(i)+0.5)*dx for i in range(nx)]
        decs  = [(float(i)+0.5)*dy-90.0 for i in range(ny)]

        # Get array geographic latitude
        geolat = self['geolat'].real()

        # Precompute latitude terms
        cos_lat = math.cos(geolat*gammalib.deg2rad)
        sin_lat = math.sin(geolat*gammalib.deg2rad)

        # Loop over all declinations and hour angles and compute the zenith
        # angle. Store the zenith angle in the map
        index = 0
        for dec in decs:
            cos_dec = math.cos(dec*gammalib.deg2rad)
            sin_dec = math.sin(dec*gammalib.deg2rad)
            for h in hours:
                cos_h  = math.cos(h*gammalib.deg2rad)
                zenith = math.acos(sin_lat*sin_dec +
                                   cos_lat*cos_dec*cos_h)*gammalib.rad2deg
                zmap[index] = zenith
                index     += 1

        # Log zenith angle map
        self._log_header2(gammalib.EXPLICIT, 'Zenith angle map')
        self._log_string(gammalib.EXPLICIT, str(zmap))

        # Return zenith angle map
        return zmap

    def _sun_radec(self, time):
        """
        Compute Right Ascension and Declination of the Sun

        Parameters
        ----------
        time : `~gammalib.GTime()`
            Time for which declination is to be computed

        Returns
        -------
        ra, dec : tuple of float
            Right Ascension and Declination of the Sun in degrees

        The formulae are inspired from
        https://en.wikipedia.org/wiki/Position_of_the_Sun
        """
        # Compute number of days since Greenwich noon, Terrestrial Time, on
        # 1 January 2000
        n = time.jd() - 2451545.0

        # Compute mean longitude of the Sun in degrees, corrected for the
        # aberration of light
        L = 280.460 + 0.9856474 * n
        while L < 0.0:
            L += 360.0
        while L >= 360.0:
            L -= 360.0

        # Compute the mean anomaly of the Sun in radians
        g = (357.528 + 0.9856003 * n) * gammalib.deg2rad

        # Compute the ecliptic longitude of the Sun in degrees
        lam = L + 1.915 * math.sin(g) + 0.020 * math.sin(2.0*g)

        # Compute Right Ascension and Declination of the Sun in degrees
        ra  = math.atan(math.cos(23.43711*gammalib.deg2rad) *
                        math.tan(lam*gammalib.deg2rad))*gammalib.rad2deg
        dec = math.asin(math.sin(23.43711*gammalib.deg2rad) *
                        math.sin(lam*gammalib.deg2rad))*gammalib.rad2deg

        # Put ra in the right quadrant (same quadrant as lam)
        while lam < 0.0:
            lam += 360.0
        while lam >= 360.0:
            lam -= 360.0
        while ra < 0.0:
            ra += 360.0
        while ra >= 360.0:
            ra -= 360.0
        q_lam = int(lam / 90.0)
        q_ra  = int(ra / 90.0)
        ra   += (q_lam-q_ra)*90.0

        # Return declination of the Sun
        return ra, dec

    def _sun_ra_exclusion(self, time):
        """
        Compute the Right Ascension difference in degrees before and after
        noon that are excluded due to the Sun constraint

        Parameters
        ----------
        time : `~gammalib.GTime()`
            Time for which the Sun exclusion is to be computed

        Returns
        -------
        dra : float
            Positive Right Ascension difference in degrees.
        """
        # Get array geographic latitude and minimum Sun zenith angle in radians
        geolat    = self['geolat'].real() * gammalib.deg2rad
        sunzenith = self['sunzenith'].real() * gammalib.deg2rad

        # Compute Right Ascension and Declination of Sun in degrees and
        # convert the Declination into radians
        sundec = self._sun_radec(time)[1] * gammalib.deg2rad

        # Compute some sines and cosines
        cos_sunzenith = math.cos(sunzenith)
        sin_geolat    = math.sin(geolat)
        cos_geolat    = math.cos(geolat)
        sin_sundec    = math.sin(sundec)
        cos_sundec    = math.cos(sundec)
        
        # Compute Right Ascension difference when Sun is below the mimimum
        # zenith angle in degrees
        dra = math.acos((cos_sunzenith - sin_geolat * sin_sundec) /
                        (cos_geolat * cos_sundec)) * gammalib.rad2deg

        # Return
        return dra

    def _hour_angle_weight(self):
        """
        Compute hour angle weights

        Computes an array specifying how many hours the array was observing
        for a given hour angle during the time interval [tmin,tmax]. The hour
        angle runs from 0 to 360 degrees.

        Compute during which time the array was observing due to hour angle
        constraints. We take here a very simple model where the night lasts
        between 6 and 10 hours (or 90-150 degrees), with a sinusoidal variation
        along the year. Every day the hour angle of the Sun set is advancing by
        4 minutes.
        """
        # Write header
        self._log_header2(gammalib.NORMAL, 'Hour angle weights')

        # Get array geographic longitude
        #geolon = self['geolon'].real()

        # Get time interval
        tmin = gammalib.GTime(self['tmin'].real(), 's')
        tmax = gammalib.GTime(self['tmax'].real(), 's')

        # Initialise hour angle list
        hour_angles = []
        
        # Set number of hour angle bins and compute conversion factor and weight
        nsteps = 1000
        ra2inx = float(nsteps)/360.0 # Conversion from RA (deg) to index
        weight = 24.0/float(nsteps)  # Weight per hour angle step
    
        # Initialise hour angle array and weight
        hours  = [0.0 for i in range(nsteps)]

        # Initialise time and loop until the end time is reached
        time = tmin
        while time <= tmax:

            # Compute Right Ascension and Declination of Sun in degrees
            sunra, sundec = self._sun_radec(time)

            # Compute by how much the zenith angle map needs to be shifted
            # to correspond to the actual time
            #
            # NOTE: The local apparent siderial time only is relevant if
            # the time interval is shorter than a day. Only in that case
            # we have to set the corresponding hour angles to zero. We need
            # to think how to properly implement that.
            #last = time.last(geolon) * 15.0
            #print(last)

            # Compute the time from now when the Sun will culminate. The time
            # is here expressed in degrees
            dra_sun = self._sun_ra_exclusion(time)
            
            # Set [0,ra_start] and [ra_stop,360.0]
            ra_start = sunra - dra_sun
            ra_stop  = sunra + dra_sun
            
            # Case 1: The RA interval during which the Sun is above the zenith
            # angle constraint is fully comprised within the [0,360] interval.
            # In that case the dark time is comprised of two intervals:
            # [0,ra_start] and [ra_stop,360]
            if ra_start >= 0.0 and ra_stop <= 360.0:
                inx_stop  = int(ra_start * ra2inx + 0.5)
                inx_start = int(ra_stop  * ra2inx + 0.5)
                for i in range(0,inx_stop):
                    hours[i] += weight
                for i in range(inx_start,nsteps):
                    hours[i] += weight

                # Log setting
                logs = 'Sun=(%.3f,%.3f) dRA=%.3f RA_excl=[%.3f,%.3f] '\
                       'Indices=[0-%d] & [%d-%d] Dark time=%.2f h %s' % \
                       (sunra, sundec, dra_sun, ra_start, ra_stop, inx_stop, inx_start,
                        nsteps, float(inx_stop+(nsteps-inx_start))*weight,
                        '(Case 1)')
                self._log_value(gammalib.VERBOSE, time.utc(), logs)

            # Case 2: The RA interval during which the Sun is above the zenith
            # angle constraint is starting at negative RA. In that case the
            # dark time is comprised of a single interval
            # [ra_stop,ra_start+360]
            elif ra_start < 0.0:
                inx_start = int(ra_stop * ra2inx + 0.5)
                inx_stop  = int((ra_start+360.0) * ra2inx + 0.5)
                for i in range(inx_start,inx_stop):
                    hours[i] += weight

                # Log setting
                logs = 'Sun=(%.3f,%.3f) dRA=%.3f RA_excl=[%.3f,%.3f] '\
                       'Indices=[%d-%d] Dark time=%.2f h %s' % \
                       (sunra, sundec, dra_sun, ra_start, ra_stop, inx_start, inx_stop,
                        float(inx_stop-inx_start)*weight, '(Case 2)')
                self._log_value(gammalib.VERBOSE, time.utc(), logs)

            # Case 3: The RA interval during which the Sun is above the zenith
            # angle constraint is stopping at RA>360. In that case the dark
            # time is comprised of a single interval [ra_stop-360,ra_start]
            elif ra_stop > 360.0:
                inx_start = int((ra_stop-360.0) * ra2inx + 0.5)
                inx_stop  = int(ra_start * ra2inx + 0.5)
                for i in range(inx_start,inx_stop):
                    hours[i] += weight

                # Log setting
                logs = 'Sun=(%.3f,%.3f) dRA=%.3f RA_excl=[%.3f,%.3f] '\
                       'Indices=[%d-%d] Dark time=%.2f h %s' % \
                       (sunra, sundec, dra_sun, ra_start, ra_stop, inx_start, inx_stop,
                        float(inx_stop-inx_start)*weight, '(Case 3)')
                self._log_value(gammalib.VERBOSE, time.utc(), logs)

            # Add seconds of one day and start with next day
            time += 86400.0

        # Build hour angle list
        dh = 360.0/float(nsteps)
        for i in range(nsteps):
            hour_angle = {'angle': dh*float(i), 'hours': hours[i]}
            hour_angles.append(hour_angle)

        # Log hour angle weights
        total = 0.0
        for hour_angle in hour_angles:
            total += hour_angle['hours']
        self._log_value(gammalib.EXPLICIT, 'Observing time', str(total)+' h')
        self._log_value(gammalib.EXPLICIT, 'Number of hour angle bins', len(hour_angles))

        # Return hour angle weights
        return hour_angles

    def _visibility_cube(self):
        """
        Compute visibility cube
        
        Compute visibility cube by displacing the zenith angle map for all
        hour angles. The visibility cube contains the number of hours a given
        celestial position is visible under a given zenith angle interval.
        Summing over all zenith angle intervals specifies for how long a given
        celestial position will be visible.
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Compute visibility cube')

        # Set visibility cube and zenith angle map dimensions and bin size
        binsz = self['binsz'].real()
        nx    = int(360.0/binsz+0.5)
        ny    = int(180.0/binsz+0.5)
        dx    = 360.0/float(nx)
        dy    = 180.0/float(ny)

        # Set zenith angle dimension and bin size
        dz   = self['dz'].real()
        zmax = self['zmax'].real()
        nz   = int(zmax/dz+0.5)
        dz   = zmax/float(nz)

        # Log dimensions
        self._log_header2(gammalib.NORMAL, 'Visibility cube dimensions')
        self._log_value(gammalib.NORMAL, 'Number of RA bins', nx)
        self._log_value(gammalib.NORMAL, 'Number of Dec bins', ny)
        self._log_value(gammalib.NORMAL, 'Number of zenith angles', nz)
        self._log_value(gammalib.NORMAL, 'RA pixel size', str(dx)+' deg')
        self._log_value(gammalib.NORMAL, 'Dec pixel size', str(dy)+' deg')
        self._log_value(gammalib.NORMAL, 'Zenith angle bin size', str(dz)+' deg')
        self._log_value(gammalib.NORMAL, 'Maximum zenith angle', str(zmax)+' deg')

        # Compute zenith angle map
        zmap = self._zenith_angle_map(nx,ny,dx,dy)

        # Setup hour angle weights. They specify for how many hours a given
        # hour angle is observed during the covered time period.
        hour_angles = self._hour_angle_weight()

        # Initialise visibility cube
        self._cube = gammalib.GSkyMap('CAR','CEL',0.0,0.0,-dx,dy,nx,ny,nz)

        # Loop over all hour angle weights
        for hour_angle in hour_angles:

            # Compute temporal shift in degrees
            shift = hour_angle['angle']

            # Initialise shifted zenith angle map
            zmap_shift = gammalib.GSkyMap('CAR','CEL',shift,0.0,-dx,dy,nx,ny)

            # Merge zenith angle map into shifted map
            zmap_shift += zmap

            # Loop over all pixels of the shifted map and add the hours during
            # which the shifted map occurs to the relevant zenith angle bin
            for i, zenith in enumerate(zmap_shift):
                iz = int(zenith/dz)
                if iz < nz:
                    self._cube[i,iz] += hour_angle['hours']

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

        # Compute visibility cube
        self._visibility_cube()

        # Optionally publish map
        if self['publish'].boolean():
            self.publish()

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
        self._log_value(gammalib.NORMAL, 'Visibility cube file', outfile.url())

        # Save the visibility cube
        self._cube.save(outfile, self['clobber'].boolean())

        # Return
        return

    def publish(self, name=''):
        """
        Publish visibility cube

        Parameters
        ----------
        name : str, optional
            Name of visibility cube
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Publish visibility cube')

        # Set default name is user name is empty
        if not name:
            user_name = self._name
        else:
            user_name = name

        # Log cube name
        self._log_value(gammalib.TERSE, 'Visibility cube name', user_name)

        # Publish cube
        self._cube.publish(user_name)

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
