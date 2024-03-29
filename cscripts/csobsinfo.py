#!/usr/bin/env python
# ==========================================================================
# Dump information about observation into log file
#
# Copyright (C) 2015-2022 Michael Mayer
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


# =============== #
# csobsinfo class #
# =============== #
class csobsinfo(ctools.csobservation):
    """
    Shows the content of an observation container
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise class members
        self._obj_dir        = None
        self._compute_offset = False
        self._offsets        = []
        self._zeniths        = []
        self._azimuths       = []
        self._pnt_ra         = []
        self._pnt_dec        = []
        self._ebounds        = gammalib.GEbounds()
        self._gti            = gammalib.GGti()

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Get parameters
        if self.obs().size() == 0:
            self._require_inobs('csobsinfo::get_parameters')
            self.obs(self._get_observations(False))

        # Initialise object position
        self._obj_dir = gammalib.GSkyDir()

        # Get (optional) offset parameters
        self._compute_offset = self['offset'].boolean()
        if self._compute_offset:
            ra  = self['ra'].real()
            dec = self['dec'].real()
            self._obj_dir.radec_deg(ra,dec)

        # Read ahead DS9 filename
        if self._read_ahead():
            self['outds9file'].query()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return


    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Initialise arrays to store certain values for reuse
        # Todo, think about using a python dictionary
        self._offsets  = []
        self._zeniths  = []
        self._azimuths = []
        self._pnt_ra   = []
        self._pnt_dec  = []
        self._ebounds  = gammalib.GEbounds()
        self._gti      = gammalib.GGti()
        obs_names      = []

        # Initialise output to be filled
        ontime         = 0.0
        livetime       = 0.0
        n_events       = 0
        n_eventbins    = 0
        n_obs_binned   = 0
        n_obs_unbinned = 0

        # Write header
        if self._logTerse():
            self._log('\n')
            self._log.header1(gammalib.number('Observation', self.obs().size()))

        # Loop over observations
        for obs in self.obs():

            # Skip non-CTA observations
            if not obs.classname() == 'GCTAObservation':
                self._log('Skipping '+obs.instrument()+' observation\n')
                continue 

            # Use observed object as observation name if name is not given
            obs_name = obs.name()
            if obs_name == '':
                obs_name = obs.object()

            # Logging
            if self._logExplicit():
                obs_id   = obs.id()
                if obs_id != '':
                    log_name = obs_name + ' (ID='+obs_id+')'
                else:
                    log_name = obs_name
                self._log.header2(log_name)

            # Retrieve observation name
            obs_names.append(obs_name)

            # Retrieve energy boundaries
            obs_bounds = obs.events().ebounds()

            # Retrieve time interval
            obs_gti = obs.events().gti()

            # Compute mean time and dead time fraction in percent
            deadfrac = (1.0-obs.deadc())*100.0

            # Retrieve pointing and store Ra,Dec
            pnt_dir = obs.pointing().dir()
            self._pnt_ra.append(pnt_dir.ra_deg())
            self._pnt_dec.append(pnt_dir.dec_deg())

            # If avaliable append energy boundaries
            if obs_bounds.size() > 0 : 
                self._ebounds.append(obs_bounds.emin(),obs_bounds.emax())

            # Append time interval
            self._gti.append(obs_gti.tstart(), obs_gti.tstop())

            # Increment global livetime and ontime
            ontime   += obs.ontime()
            livetime += obs.livetime()

            # Bookkeeping
            if obs.eventtype() == 'CountsCube':
                n_eventbins  += obs.events().size()
                n_obs_binned += 1
                is_binned     = 'yes'
                is_what       = 'Number of bins'
            else:
                n_events       +=  obs.events().size()
                n_obs_unbinned += 1
                is_binned     = 'no'
                is_what       = 'Number of events'
            self._log_value(gammalib.EXPLICIT, 'Binned', is_binned)
            self._log_value(gammalib.EXPLICIT, is_what, obs.events().size())

            # Retrieve zenith and azimuth and store for later use
            zenith  = obs.pointing().zenith()
            azimuth = obs.pointing().azimuth()
            self._zeniths.append(zenith)
            self._azimuths.append(azimuth)

            # Optionally compute offset with respect to target direction
            if self._compute_offset:
                offset = pnt_dir.dist_deg(self._obj_dir)
                self._offsets.append(offset)

            # Optionally log details
            if self._logExplicit():

                # Log the observation energy range (if available)
                self._log.parformat('Energy range')
                if obs_bounds.size() == 0:
                    self._log('undefined')
                else:
                    self._log(str(obs_bounds.emin()))
                    self._log(' - ')
                    self._log(str(obs_bounds.emax()))
                self._log('\n')

                # Log observation time interval
                self._log.parformat('Time range (MJD)')
                if obs_gti.size() == 0:
                    self._log('undefined')
                else:
                    self._log(str(obs_gti.tstart().mjd()))
                    self._log(' - ')
                    self._log(str(obs_gti.tstop().mjd()))
                self._log('\n')

            # Log observation information
            self._log_value(gammalib.EXPLICIT, 'Ontime', '%.3f s' %
                            obs.ontime())
            self._log_value(gammalib.EXPLICIT, 'Livetime', '%.3f s' %
                            obs.livetime())
            self._log_value(gammalib.EXPLICIT, 'Deadtime fraction', '%.3f %%' %
                            deadfrac)
            self._log_value(gammalib.EXPLICIT, 'Pointing', pnt_dir)

            # Optionally log offset with respect to target direction
            if self._compute_offset:
                self._log_value(gammalib.EXPLICIT,
                                'Offset from target', '%.2f deg' % offset)

            # Log Zenith and Azimuth angles
            self._log_value(gammalib.EXPLICIT, 'Zenith angle', '%.2f deg' %
                            zenith)
            self._log_value(gammalib.EXPLICIT, 'Azimuth angle', '%.2f deg' %
                            azimuth)

        # Write summary header
        self._log_header1(gammalib.NORMAL, 'Summary')

        # Log general summary
        self._log_header3(gammalib.NORMAL, 'Observations')
        self._log_value(gammalib.NORMAL, 'Unbinned observations', n_obs_unbinned)
        self._log_value(gammalib.NORMAL, 'Binned observations', n_obs_binned)
        self._log_header3(gammalib.NORMAL, 'Events')
        self._log_value(gammalib.NORMAL, 'Number of events', n_events)
        self._log_value(gammalib.NORMAL, 'Number of bins', n_eventbins)

        # Compute mean offset, azimuth and zenith angle
        if len(self._offsets) > 0:
            mean_offset = '%.2f deg' % (sum(self._offsets) / len(self._offsets))
        else:
            mean_offset = 'Unknown'
        if len(self._zeniths) > 0:
            mean_zenith = '%.2f deg' % (sum(self._zeniths) / len(self._zeniths))
        else:
            mean_zenith = 'Unknown'
        if len(self._azimuths) > 0:
            mean_azimuth = '%.2f deg' % (sum(self._azimuths) / len(self._azimuths))
        else:
            mean_azimuth = 'Unknown'

        # Log mean offset, azimuth and zenith angle
        self._log_header3(gammalib.NORMAL, 'Pointings')
        self._log_value(gammalib.NORMAL, 'Mean offset angle',  mean_offset)
        self._log_value(gammalib.NORMAL, 'Mean zenith angle',  mean_zenith)
        self._log_value(gammalib.NORMAL, 'Mean azimuth angle', mean_azimuth)

        # Optionally log names of observations. Note that the set class is
        # used to extract all different observation names from the list of
        # observation names, and the set class is only available from
        # Python 2.4 on.
        if sys.version_info >= (2,4):
            obs_set = set(obs_names)
            for name in obs_set:
                self._log_value(gammalib.EXPLICIT,'"'+name+'"',
                                obs_names.count(name))

        # Get energy boundary information
        if self._ebounds.size() == 0:
            min_value = 'undefined'
            max_value = 'undefined'
        else:
            min_value = str(self._ebounds.emin())
            max_value = str(self._ebounds.emax())

        # Log energy range
        self._log_header3(gammalib.NORMAL, 'Energy range')
        self._log_value(gammalib.NORMAL, 'Minimum energy', min_value)
        self._log_value(gammalib.NORMAL, 'Maximum energy', max_value)

        # Log time range
        mjd = '%.3f - %.3f' % (self._gti.tstart().mjd(),  self._gti.tstop().mjd())
        utc = '%s - %s'     % (self._gti.tstart().utc(),  self._gti.tstop().utc())
        self._log_header3(gammalib.NORMAL, 'Time range')
        self._log_value(gammalib.NORMAL, 'MJD (days)', mjd)
        self._log_value(gammalib.NORMAL, 'UTC',  utc)

        # Log ontime and livetime in different units
        on_time   = '%.2f s = %.2f min = %.2f h' % \
                    (ontime, ontime/60., ontime/3600.)
        live_time = '%.2f s = %.2f min = %.2f h' % \
                    (livetime, livetime/60., livetime/3600.)
        self._log_value(gammalib.NORMAL, 'Total ontime', on_time)
        self._log_value(gammalib.NORMAL, 'Total livetime', live_time)

        # Return
        return

    def save(self):
        """ 
        Save pointings into DS9 region file

        This method saves all pointing directions that are found in the
        observation container into a DS9 region file. If "NONE" is
        specified for the "outds9file" parameter the method does nothing.
        """
        # Check if DS9 file is valid
        if self['outds9file'].is_valid():

            # Get output filename in case it was not read ahead
            ds9file = self['outds9file'].filename()

            # Write header
            self._log_header1(gammalib.TERSE, 'Save pointings in DS9 file')

            # Log filename
            self._log_value(gammalib.NORMAL, 'DS9 filename', ds9file.url())

            # Open file
            f = open(ds9file.url(),'w')

            # Write coordinate system
            f.write('fk5\n')

            # Loop over pointings
            for i in range(len(self._pnt_ra)):

                # Create string
                line  = 'point('
                line += str(self._pnt_ra[i])+','+str(self._pnt_dec[i])+')'
                line += ' # point=cross 20 width=3\n'

                # Write to file
                f.write(line)

            # Close file
            f.close()

        # Return
        return

    def zeniths(self):
        """
        Return zenith angles
        """
        return self._zeniths

    def azimuths(self):
        """
        Return azimuth angles
        """
        return self._azimuths

    def ras(self):
        """
        Return pointings right ascension
        """
        return self._pnt_ra

    def decs(self):
        """
        Return pointings declination
        """
        return self._pnt_dec

    def offsets(self):
        """
        Return offset angles
        """
        return self._offsets

    def ebounds(self):
        """
        Return energy boundaries
        """
        return self._ebounds

    def gti(self):
        """
        Return good time intervals
        """
        return self._gti


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csobsinfo(sys.argv)

    # Execute application
    app.execute()
