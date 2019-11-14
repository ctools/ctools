#! /usr/bin/env python
# ==========================================================================
# Compute a visibility cube
#
# Copyright (C) 2016-2019 Juergen Knoedlseder
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
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._cube    = gammalib.GSkyMap()
        self._results = []

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get all parameters
        """
        # Query parameters
        self['mjdref'].real()
        self['tmin'].time()
        self['tmax'].time()
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

        The zenith angle of a position (ra,dec) depends on the declination and
        the hour angle h and is given by

        zenith(h,dec) = arccos( sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(h) )

        The hour angle h (or local hour angle, LHA) is defined as the difference
        between local siderial time (LST) and the Right Ascension

        h = LST - ra

        The map is computed for h=-ra which is equivalent to GST=lon (or LST=0).
        In other words, the map corresponds to the time when ra=0 passes through
        the local meridian.

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
        """
        # Initialise zenith angle map
        zmap = gammalib.GSkyMap('CAR','CEL',0.0,0.0,-dx,dy,nx,ny)

        # Set hour angle and declination vectors
        hours = [(float(i)+0.5)*dx      for i in range(nx)]
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
                index      += 1

        # Log zenith angle map
        self._log_header2(gammalib.EXPLICIT, 'Zenith angle map')
        self._log_string(gammalib.EXPLICIT, str(zmap))

        # Return zenith angle map
        return zmap

    def _sun_ra_exclusion(self, time):
        """
        Compute the half length of the Right Ascension interval to exclude due
        to the Sun constraint

        The Sun zenith angle constraint implies an interval in Right Ascension
        that is to be excluded (that's the interval during which it is day).
        This method computes half the length of this interval in degrees.

        Parameters
        ----------
        time : `~gammalib.GTime()`
            Time for which the Sun exclusion is to be computed

        Returns
        -------
        dra : float
            Right Ascension interval half length (degrees)
        """
        # Get array geographic latitude and minimum Sun zenith angle in radians
        geolat    = self['geolat'].real()    * gammalib.deg2rad
        sunzenith = self['sunzenith'].real() * gammalib.deg2rad

        # Compute Right Ascension and Declination of Sun in degrees and
        # convert the Declination into radians
        sun = gammalib.GSkyDir()
        sun.sun(time)
        sundec = sun.dec()

        # Compute some sines and cosines
        cos_sunzenith = math.cos(sunzenith)
        sin_geolat    = math.sin(geolat)
        cos_geolat    = math.cos(geolat)
        sin_sundec    = math.sin(sundec)
        cos_sundec    = math.cos(sundec)

        # Compute Right Ascension difference when Sun is below the mimimum
        # zenith angle in degrees
        dra   = 0.0
        nom   = cos_sunzenith - sin_geolat * sin_sundec
        denom = cos_geolat * cos_sundec
        if denom != 0.0:
            arg = nom/denom
            if arg >= -1.0 and arg <= 1.0:
                dra = math.acos(arg) * gammalib.rad2deg

        # Return
        return dra

    def _moon_ra_exclusion(self, time):
        """
        Compute the half length of the Right Ascension interval to exclude due
        to the Moon constraint
        
        The Moon zenith angle constraint implies an interval in Right Ascension
        that is to be excluded. This method computes half the length of this
        interval.

        Parameters
        ----------
        time : `~gammalib.GTime()`
            Time for which the Moon exclusion is to be computed

        Returns
        -------
        dra : float
            Right Ascension interval half length (degrees)
        """
        # Get array geographic latitude and minimum Moon zenith angle in radians
        geolat     = self['geolat'].real()     * gammalib.deg2rad
        moonzenith = self['moonzenith'].real() * gammalib.deg2rad

        # Compute Right Ascension and Declination of Moon in degrees and
        # convert the Declination into radians
        moon = gammalib.GSkyDir()
        moon.moon(time)
        moondec = moon.dec()

        # Compute some sines and cosines
        cos_moonzenith = math.cos(moonzenith)
        sin_geolat     = math.sin(geolat)
        cos_geolat     = math.cos(geolat)
        sin_moondec    = math.sin(moondec)
        cos_moondec    = math.cos(moondec)

        # Compute Right Ascension difference when Moon is below the mimimum
        # zenith angle in degrees
        dra   = 0.0
        nom   = cos_moonzenith - sin_geolat * sin_moondec
        denom = cos_geolat * cos_moondec
        if denom != 0.0:
            arg = nom/denom
            if arg >= -1.0 and arg <= 1.0:
                dra = math.acos(arg) * gammalib.rad2deg

        # Return
        return dra

    def _adjust_ra_interval(self, ra_start, ra_stop):
        """
        Adjust Right Ascension interval so that it overlaps with the [0,360[
        interval

        Parameters
        ----------
        ra_start : float
            Start of interval to exclude (degrees)
        ra_stop : float
            Stop of interval to exclude (degrees)

        Returns
        -------
        ra_start : float
            Adjusted start of interval to exclude (degrees)
        ra_stop : float
            Adjusted stop of interval to exclude (degrees)
        """
        # Adjust interval if both boundaries are negative
        while ra_start < 0.0 and ra_stop < 0.0:
            ra_start += 360.0
            ra_stop  += 360.0

        # Adjust interval if both boundaries are equal or larger than 360
        # degrees
        while ra_start >= 360.0 and ra_stop >= 360.0:
            ra_start -= 360.0
            ra_stop  -= 360.0

        # Return interval
        return ra_start, ra_stop

    def _exclude_ra_interval(self, hours, ra_start, ra_stop):
        """
        Exclude Right Ascension interval from array of hour angles

        Parameters
        ----------
        hours : list of floats
            Array of hours
        ra_start : float
            Start of interval to exclude (degrees)
        ra_stop : float
            Stop of interval to exclude (degrees)

        Returns
        -------
        hours : list of floats
            Array of hours with excluded Right Ascension interval
        """
        # Conversion from RA (degrees) to index
        nhours = len(hours)
        ra2inx = float(nhours)/360.0

        # Adjust interval
        ra_start, ra_stop = self._adjust_ra_interval(ra_start, ra_stop)

        # Case 1: The RA interval is fully comprised within the [0,360]
        #         interval, hence we simply exclude this interval.
        if ra_start >= 0.0 and ra_stop <= 360.0:
            inx_start = int(ra_start * ra2inx + 0.5)
            inx_stop  = int(ra_stop  * ra2inx + 0.5)
            for i in range(inx_start,inx_stop):
                hours[i] = 0.0

        # Case 2: The RA interval is starting at negative RA. In that case
        #         there are two intervals to exclude:
        #         [ra_start+360,360] and [0,ra_stop]
        elif ra_start < 0.0:
            inx_start = int((ra_start+360.0) * ra2inx + 0.5)
            inx_stop  = int(ra_stop * ra2inx + 0.5)
            for i in range(inx_start,nhours):
                hours[i] = 0.0
            for i in range(0,inx_stop):
                hours[i] = 0.0

        # Case 3: The RA interval is stopping at RA>360. In that case there
        #         are two intervals to exclude:
        #         [0,ra_stop-360] and [0,ra_start]
        elif ra_stop > 360.0:
            inx_start = int(ra_start * ra2inx + 0.5)
            inx_stop  = int((ra_stop-360.0) * ra2inx + 0.5)
            for i in range(0,inx_stop):
                hours[i] = 0.0
            for i in range(inx_start,nhours):
                hours[i] = 0.0

        # Return hours
        return hours

    def _hour_angle_weight(self):
        """
        Compute hour angle weights

        Computes an array that specifies the time during which a given hour
        angle was observed during the observing time interval [tmin,tmax].
        The hour angle runs from 0 to 360 degrees, array values are in units
        of hours.

        The method loops over the days in the time interval [tmin,tmax], with
        time = tmin + i*86400 seconds. At each time step the position in Right
        Ascension and Declination of the Sun is computed (the Sun's
        Declination is in fact not used here).

        An interval of [ra_sun-dra_sun, ra_sun+dra_sun] is assumed to be the
        day. The method _sun_ra_exclusion() is used to compute dra_sun, which
        is half of the length of the day in degrees (recall that 15 degrees
        is one hour). The length of the day is computed using the sunzenith
        constraint.
        """
        # Write header
        self._log_header2(gammalib.NORMAL, 'Hour angle weights')

        # Get MET time reference
        tref = gammalib.GTimeReference(self['mjdref'].real(),'s','TT','LOCAL')

        # Get time interval
        tmin = self['tmin'].time(tref)
        tmax = self['tmax'].time(tref)

        # Initialise hour angle list and results
        hour_angles   = []
        self._results = []

        # Set number of hour angle bins and compute conversion factor and weight
        nsteps = 1000
        ra2inx = float(nsteps)/360.0 # Conversion from RA (deg) to index
        weight = 24.0/float(nsteps)  # Weight per hour angle step

        # Initialise hour angle array and weight
        hours  = [0.0 for i in range(nsteps)]

        # Initialise time and loop until the end time is reached
        time = tmin
        while time <= tmax:

            # Initialise hours array for this time step
            hours_time = [weight for i in range(nsteps)]

            # Compute half the length of the exclusion intervals in degrees
            sun_dra  = self._sun_ra_exclusion(time)
            moon_dra = self._moon_ra_exclusion(time)

            # Compute Right Ascension and Declination of Sun in degrees and
            # derive Sun exclusion interval
            sun = gammalib.GSkyDir()
            sun.sun(time)
            sun_ra       = sun.ra_deg()
            sun_dec      = sun.dec_deg()
            sun_ra_start = sun_ra - sun_dra
            sun_ra_stop  = sun_ra + sun_dra

            # Adjust interval
            sun_ra_start, sun_ra_stop = self._adjust_ra_interval(sun_ra_start,
                                                                 sun_ra_stop)

            # Compute Right Ascension and Declination of Moon in degrees and
            # derive Moon exclusion interval
            moon = gammalib.GSkyDir()
            moon.moon(time)
            moon_ra       = moon.ra_deg()
            moon_dec      = moon.dec_deg()
            moon_ra_start = moon_ra - moon_dra
            moon_ra_stop  = moon_ra + moon_dra

            # Adjust interval
            moon_ra_start, moon_ra_stop = self._adjust_ra_interval(moon_ra_start,
                                                                   moon_ra_stop)

            # Compute Moon elongation and illumation fraction (Moon phase)
            elongation = sun.dist_deg(moon)
            fli        = (1.0 - math.cos(elongation * gammalib.deg2rad))/2.0

            # Exclude hours due to Sun constraint (this defines the night)
            hours_time = self._exclude_ra_interval(hours_time, sun_ra_start,
                                                               sun_ra_stop)

            # Exclude hours due to Moon constraint if the illumination fraction
            # is equal to or above the maximim fraction of illumination
            if fli >= self['maxfli'].real():
                hours_time = self._exclude_ra_interval(hours_time, moon_ra_start,
                                                                   moon_ra_stop)

            # Add hours
            dark_time = 0.0
            for i in range(nsteps):
                hours[i]  += hours_time[i]
                dark_time += hours_time[i]

            # Set result record
            result = {'time': time.copy(),
                      'sun_ra': sun_ra, 'sun_dec': sun_dec,
                      'moon_ra': moon_ra, 'moon_dec': moon_dec,
                      'elongation': elongation, 'fli': fli,
                      'sun_ra_start': sun_ra_start, 'sun_ra_stop': sun_ra_stop,
                      'moon_ra_start': moon_ra_start, 'moon_ra_stop': moon_ra_stop,
                      'dark_time' : dark_time}

            # Append result record to results
            self._results.append(result)

            # Log results
            logs = 'Sun=(%8.3f,%7.3f) Moon=(%8.3f,%7.3f) '\
                   'Sun_RA_excl=[%8.3f,%8.3f] Moon_RA_excl=[%8.3f,%8.3f] '\
                   'FLI=%4.2f Dark time=%5.2f h' % \
                   (sun_ra, sun_dec, moon_ra, moon_dec,
                    sun_ra_start, sun_ra_stop, moon_ra_start, moon_ra_stop,
                    fli, dark_time)
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

    def _save_results(self, outfile, clobber):
        """
        Save results in VISIBILITY FITS table

        Parameters
        ----------
        outfile : str
            Result FITS file name
        clobber : bool
            Overwrite existing file?
        """
        # Create FITS table columns
        nrows         = len(self._results)
        time          = gammalib.GFitsTableStringCol('Time', nrows, 20)
        mjd           = gammalib.GFitsTableDoubleCol('MJD', nrows)
        sun_ra        = gammalib.GFitsTableDoubleCol('RA_sun', nrows)
        sun_dec       = gammalib.GFitsTableDoubleCol('DEC_sun', nrows)
        moon_ra       = gammalib.GFitsTableDoubleCol('RA_moon', nrows)
        moon_dec      = gammalib.GFitsTableDoubleCol('DEC_moon', nrows)
        sun_ra_start  = gammalib.GFitsTableDoubleCol('RA_sun_start', nrows)
        sun_ra_stop   = gammalib.GFitsTableDoubleCol('RA_sun_stop', nrows)
        moon_ra_start = gammalib.GFitsTableDoubleCol('RA_moon_start', nrows)
        moon_ra_stop  = gammalib.GFitsTableDoubleCol('RA_moon_stop', nrows)
        elongation    = gammalib.GFitsTableDoubleCol('Elongation', nrows)
        fli           = gammalib.GFitsTableDoubleCol('FLI', nrows)
        dark_time     = gammalib.GFitsTableDoubleCol('Darktime', nrows)

        # Set units of table columns
        mjd.unit('days')
        sun_ra.unit('deg')
        sun_dec.unit('deg')
        moon_ra.unit('deg')
        moon_dec.unit('deg')
        sun_ra_start.unit('deg')
        sun_ra_stop.unit('deg')
        moon_ra_start.unit('deg')
        moon_ra_stop.unit('deg')
        elongation.unit('deg')
        dark_time.unit('hours')

        # File FITS table columns
        for i, result in enumerate(self._results):
            time[i]          = result['time'].utc()
            mjd[i]           = result['time'].mjd()
            sun_ra[i]        = result['sun_ra']
            sun_dec[i]       = result['sun_dec']
            moon_ra[i]       = result['moon_ra']
            moon_dec[i]      = result['moon_dec']
            sun_ra_start[i]  = result['sun_ra_start']
            sun_ra_stop[i]   = result['sun_ra_stop']
            moon_ra_start[i] = result['moon_ra_start']
            moon_ra_stop[i]  = result['moon_ra_stop']
            elongation[i]    = result['elongation']
            fli[i]           = result['fli']
            dark_time[i]     = result['dark_time']

        # Initialise FITS Table with extension "VISIBILITY"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('VISIBILITY')

        # Add Header for compatibility with gammalib.GMWLSpectrum
        table.card('INSTRUME', 'CTA', 'Name of Instrument')
        table.card('TELESCOP', 'CTA', 'Name of Telescope')

        # Append filled columns to fits table    
        table.append(time)
        table.append(mjd)
        table.append(sun_ra)
        table.append(sun_dec)
        table.append(moon_ra)
        table.append(moon_dec)
        table.append(sun_ra_start)
        table.append(sun_ra_stop)
        table.append(moon_ra_start)
        table.append(moon_ra_stop)
        table.append(elongation)
        table.append(fli)
        table.append(dark_time)

        # Append table to result FITS file
        fits = gammalib.GFits(outfile)
        fits.append(table)
        fits.save(clobber)

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

        # Save the results in VISIBILITY table extension
        self._save_results(outfile, self['clobber'].boolean())

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
            user_name = self._name()
        else:
            user_name = name

        # Log cube name
        self._log_value(gammalib.TERSE, 'Visibility cube name', user_name)

        # Publish cube
        self._cube.publish(user_name)

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
