#! /usr/bin/env python
# ==========================================================================
# Generate COMPTEL database from HEASARC archive
#
# Copyright (C) 2022-2023 Juergen Knoedlseder
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
from distutils.spawn import find_executable
import gammalib
import ctools
from cscripts import modutils


# ============== #
# comgendb class #
# ============== #
class comgendb(ctools.cscript):
    """
    Generate COMPTEL database from HEASARC archive
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._evp_list = []
        self._oad_list = []
        self._hkd_list = []
        self._tim_list = []
        self._xml_list = []

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query parameters
        self['archive'].string()
        self['refdata'].query()
        self['quality'].integer()
        self['download'].boolean()

        # Query ahead output model filename
        if self._read_ahead():
            self['dbase'].string()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _download(self):
        """
        Download COMPTEL HEASARC archive
        """
        # Log header
        self._log_header1(gammalib.TERSE, 'Download HEASARC archive')

        # Check whether wget executable exists
        if find_executable('wget') is None:
            msg = ('"wget" command not found on system, unable to download '
                   'COMPTEL HEASARC archive.')
            raise RuntimeError(msg)

        # Set archive directory name
        dirname = gammalib.expand_env(self['archive'].string())

        # Build wget command string
        cmd = ('wget -r -N -nH --cut-dirs=4 -A "*.fits,*.fits.gz" -e robots=off '
               '-P %s -q https://heasarc.gsfc.nasa.gov/FTP/compton/data/comptel/' %
               (dirname))

        # Launch download
        os.system(cmd)

        # Return
        return

    def _get_phases(self):
        """
        Get available phases
        """
        # Set archive directory name
        dirname = gammalib.expand_env(self['archive'].string())

        # Raise an exception if the archive does not exist
        if not os.path.isdir(dirname):
            msg = ('Archive directory "%s" does not exist.' % dirname)
            raise RuntimeError(msg)

        # Set search pattern
        pattern = '%s/phase*' % (gammalib.expand_env(self['archive'].string()))

        # Get available phases
        phases = glob.glob(pattern)
        phases.sort()

        # Return phases
        return phases

    def _get_relpath(self, path, dest):
        """
        Get relative path with respect to dest

        Parameters
        ----------
        path : str
            Full path name
        dest : str
            Destination path name

        Returns
        -------
        relpath : str
            Relative path
        """
        # Split path and destination in segments
        path_split = path.split('/')
        dest_split = dest.split('/')

        # Get index of last agreement
        index = -1
        for i, segment in enumerate(dest_split):
            if segment == path_split[i]:
                index = i

        # Set relative path
        if index == -1:
            relpath = path
        else:
            relpath = ''
            for i in range(index+1,len(path_split)):
                if relpath != '':
                    relpath += '/'
                relpath += path_split[i]

        # Return relative path
        return relpath

    def _sorted(self, elements):
        """
        Sort list of dict according to 'tjd_start'

        Parameters
        ----------
        elements : list of dict
            List of dictionary

        Returns
        -------
        sorted_elements : list of dict
            Sorted list of dictionary
        """
        # Set Python version that has sorted function and current
        # Python version
        req_version = (2,4)
        cur_version = sys.version_info

        # If Python version has sorted the use sorted function
        if cur_version >= req_version:
            sorted_elements = sorted(elements, key=lambda k: k['tjd_start'])

        # ... otherwise use sort() with cmp
        else:
            sorted_elements = elements
            sorted_elements.sort(lambda x,y : cmp(x['tjd_start'], y['tjd_start']))

        # Return sorted elements
        return sorted_elements

    def _get_evp_list(self):
        """
        Get EVP list

        Returns
        -------
        evp_list : list of dict
            EVP database
        """
        # Initialse EVP list
        evp_list = []

        # Get available phases
        phases = self._get_phases()

        # Loop over phases
        for phase in phases:

            # Get available VPs
            vps = glob.glob('%s/vp*' % phase)
            vps.sort()

            # Loop over VPs
            for vp in vps:

                # Get EVPs for VP
                evps = self._get_evps_for_vp(vp)

                # Append EVPs to list
                evp_list.extend(evps)

        # Sort EVP list by start day
        evp_list = self._sorted(evp_list)

        # Return EVP list
        return evp_list

    def _get_evps_for_vp(self, vp):
        """
        Build EVP database entries for Viewing Period

        Parameters
        ----------
        vp : str
            VP path

        Returns
        -------
        evp_list : list of dict
            EVP database entries
        """
        # Get archive name and expanded path
        archive      = self['archive'].string()
        archive_path = gammalib.expand_env(archive)

        # Get VP name
        vpname = os.path.basename(vp)

        # Initialse EVP list
        evp_list = []

        # Get EVP files
        evps = glob.glob('%s/*evp.fits*' % vp)

        # Loop over all EVP files
        for evp in evps:

            # Set EVP filename
            evpname = archive+'/'+self._get_relpath(evp, archive_path).strip('.gz')

            # Try opening file
            try:
                fits = gammalib.GFits(evpname)
            except:
                self._log_string(gammalib.NORMAL, '*** Unable to open EVP file'
                                 ' "%s". Skip file.' % (evpname))
                continue

            # Get information from EVP file
            evp_dict = self._get_evp_info(evpname)

            # Skip file if it has a negative quality
            if evp_dict['quality'] < 0:
                self._log_string(gammalib.NORMAL, '*** EVP file "%s" has negative '
                                 'quality %d. Skip file.' % (evpname, evp_dict['quality']))
                continue

            # Skip file if the duration derived from the events is not consistent
            # with the duration derived from the header
            duration_evt = evp_dict['tjd_stop'] - evp_dict['tjd_start']
            duration_hdr = evp_dict['vieday']   - evp_dict['visday']
            difference   = abs(duration_evt-duration_hdr)
            if difference > 3:
                self._log_string(gammalib.NORMAL, '*** EVP file "%s" has inconsistent '
                                 'duration between events (%d days) and header '
                                 '(%d days). Skip file.' % \
                                 (evpname, duration_evt, duration_hdr))
                continue
            if duration_hdr > duration_evt:
                self._log_string(gammalib.NORMAL, '--- EVP file "%s" header indicates '
                                 'duration of %d days, but events cover only %d days.' % \
                                 (evpname, duration_hdr, duration_evt))
            elif duration_hdr < duration_evt:
                self._log_string(gammalib.NORMAL, '--- EVP file "%s" header indicates '
                                 'duration of %d days, but events cover %d days.' % \
                                 (evpname, duration_hdr, duration_evt))

            # Set VP
            evp_dict['vp'] = vpname

            # Append dictionary
            evp_list.append(evp_dict)

        # Sort EVP list by start day
        evp_list = self._sorted(evp_list)

        # Return EVP list
        return evp_list

    def _get_evp_info(self, evpname):
        """
        Build EVP database entry from EVP file

        Parameters
        ----------
        evpname : str
            EVP filename

        Returns
        -------
        evp_dict : dict
            EVP database entry
        """
        # Open EVP file
        fits = gammalib.GFits(evpname)

        # Get event table
        events  = fits.table(1)
        nevents = events.nrows()

        # Extract attributes for EVP file
        obs_id   = events.real('OBS_ID')
        object   = events.string('OBJECT')
        visday   = events.integer('VISDAY')
        vistim   = events.integer('VISTIM')
        vieday   = events.integer('VIEDAY')
        vietim   = events.integer('VIETIM')
        glon_scz = events.real('GLON_SCZ')
        glat_scz = events.real('GLAT_SCZ')
        verno    = events.integer('DSD_REP')
        quality  = events.integer('DSD_QUA')
        if events.has_card('TITLE'):
            title = events.string('TITLE')
        elif events.has_card('TITLE1') and events.has_card('TITLE2'):
            title = events.string('TITLE1') + events.string('TITLE2')
        dsd_id = '%s-%s-%d' % (events.string('DSD_SB'), events.string('DSTID'),
                               events.integer('DSD_NO'))

        # Get time columns
        tjd  = events['TJD']
        tics = events['TICS']

        # Get start and end time. Catch any errors since files may be corrupt.
        # In case of errors, use times from file header.
        try:
            tjd_start  = tjd[0]
            tics_start = tics[0]
            tjd_stop   = tjd[nevents-1]
            tics_stop  = tics[nevents-1]
        except:
            self._log_string(gammalib.NORMAL, '*** Corrupt EVP file "%s" encountered. '
                             'Set EVP quality to -200.' % (evpname))
            tjd_start  = visday
            tics_start = vistim
            tjd_stop   = vieday
            tics_stop  = vietim
            quality    = -200

        # Compute GTime of first and last event
        tstart = gammalib.com_time(tjd_start, tics_start)
        tstop  = gammalib.com_time(tjd_stop, tics_stop)

        # Check that start time is before stop time
        if tstart >= tstop:
            quality = -200
            gti     = gammalib.GGti()
            self._log_string(gammalib.NORMAL, '*** Corrupt EVP file "%s" encountered. '
                             'Start time %s is later than stop time %s. Set EVP '
                             'quality to -200.' % (evpname, tstart.utc(), tstop.utc()))

        # Set GTI for EVP
        else:
            gti = gammalib.GGti(tstart, tstop)

        # Build dictionary
        evp_dict = {'filename':   evpname,
                    'nevents':    nevents,
                    'obs_id':     obs_id,
                    'dsd_id':     dsd_id,
                    'object':     object,
                    'glon_scz':   glon_scz,
                    'glat_scz':   glat_scz,
                    'date_start': tstart.utc(),
                    'date_stop':  tstop.utc(),
                    'tjd_start':  tjd_start,
                    'tics_start': tics_start,
                    'tjd_stop':   tjd_stop,
                    'tics_stop':  tics_stop,
                    'tstart':     tstart,
                    'tstop':      tstop,
                    'visday':     visday,
                    'vistim':     vistim,
                    'vieday':     vieday,
                    'vietim':     vietim,
                    'gti':        gti,
                    'quality':    quality,
                    'verno':      verno,
                    'title':      title}

        # Close FITS file
        fits.close()

        # Correct quality flags for some files according to database reports
        if evp_dict['dsd_id'] == 'M-EVP-26342' or \
           evp_dict['dsd_id'] == 'M-EVP-24955' or \
           evp_dict['dsd_id'] == 'M-EVP-24963':
            self._log_string(gammalib.NORMAL, '*** Quality flag for "%s" changed '
                             'from %s to 50 according to database reports.' %
                             (evp_dict['filename'], evp_dict['quality']))
            evp_dict['quality'] = 50

        # Return information
        return evp_dict

    def _get_oad_list(self):
        """
        Get OAD list

        Returns
        -------
        oad_list : list of dict
            OAD database
        """
        # Initialse OAD list
        oad_list = []

        # Get available phases
        phases = self._get_phases()

        # Loop over phases
        for phase in phases:

            # Get available VPs
            vps = glob.glob('%s/vp*' % phase)
            vps.sort()

            # Loop over VPs
            for vp in vps:

                # Get OADs for VP
                oads = self._get_oads_for_vp(vp)

                # Append OADs to list
                oad_list.extend(oads)

        # Sort OAD list by start day
        oad_list = self._sorted(oad_list)

        # Return OAD list
        return oad_list

    def _get_oads_for_vp(self, vp):
        """
        Build OAD database entries for Viewing Period

        Parameters
        ----------
        vp : str
            VP path

        Returns
        -------
        oad_list : list of dict
            OAD database entries
        """
        # Get archive name and expanded path
        archive      = self['archive'].string()
        archive_path = gammalib.expand_env(archive)

        # Get VP name
        vpname = os.path.basename(vp)

        # Initialse OAD list
        oad_list = []

        # Get OAD files
        oads = glob.glob('%s/*oad.fits*' % vp)

        # Loop over all OAD files
        for oad in oads:

            # Set OAD filename
            oadname = archive+'/'+self._get_relpath(oad, archive_path).strip('.gz')

            # Try opening file
            try:
                fits = gammalib.GFits(oadname)
            except:
                self._log_string(gammalib.NORMAL, '*** Unable to open OAD file "%s". '
                                 'Skip file.' % (oadname))
                continue

            # Get information from OAD file
            oad_dict = self._get_oad_info(oadname)

            # Skip file if it has a negative quality
            if oad_dict['quality'] < 0:
                self._log_string(gammalib.NORMAL, '*** OAD file "%s" has quality %d. '
                                 'Skip file.' % (oadname, oad_dict['quality']))
                continue

            # Set VP
            oad_dict['vp'] = vpname

            # Append dictionary
            oad_list.append(oad_dict)

        # Sort OAD list by start day
        oad_list = self._sorted(oad_list)

        # Return OAD list
        return oad_list

    def _get_oad_info(self, oadname):
        """
        Build OAD database entry from OAD file

        Parameters
        ----------
        oadname : str
            OAD filename

        Returns
        -------
        oad_dict : dict
            OAD database entry
        """
        # Open OAD file
        fits = gammalib.GFits(oadname)

        # Get superpacket table
        records  = fits.table(1)
        nrecords = records.nrows()

        # Extract attributes for OAD file
        obs_id   = records.real('OBS_ID')
        object   = records.string('OBJECT')
        if records.has_card('GLON_SCZ'):
            glon_scz = records.real('GLON_SCZ')
        else:
            glon_scz = -1000.0
        if records.has_card('GLAT_SCZ'):
            glat_scz = records.real('GLAT_SCZ')
        else:
            glat_scz = -1000.0
        quality = records.integer('DSD_QUA')
        dsd_id  = '%s-%s-%d' % (records.string('DSD_SB'), records.string('DSTID'),
                                records.integer('DSD_NO'))

        # Get time columns
        tjd  = records['TJD']
        tics = records['TICS']

        # Get start and end time. Catch any errors since files may be corrupt.
        # In case of errors, use times from file header.
        try:
            tjd_start  = tjd[0]
            tics_start = tics[0]
            tjd_stop   = tjd[nrecords-1]
            tics_stop  = tics[nrecords-1]
        except:
            self._log_string(gammalib.NORMAL, '*** Corrupt OAD file "%s" encountered. '
                             'Set OAD quality to -200.' % (oadname))
            tjd_start  = records.integer('VISDAY')
            tics_start = records.integer('VISTIM')
            tjd_stop   = records.integer('VIEDAY')
            tics_stop  = records.integer('VIETIM')
            quality    = -200

        # Compute GTime of first and last superpacket
        tstart = gammalib.com_time(tjd_start, tics_start)
        tstop  = gammalib.com_time(tjd_stop, tics_stop)

        # Set GTI for OAD
        gti = gammalib.GGti(tstart, tstop)

        # Build dictionary
        oad_dict = {'filename':     oadname,
                    'superpackets': nrecords,
                    'obs_id':       obs_id,
                    'dsd_id':       dsd_id,
                    'object':       object,
                    'glon_scz':     glon_scz,
                    'glat_scz':     glat_scz,
                    'date_start':   tstart.utc(),
                    'date_stop':    tstop.utc(),
                    'tjd_start':    tjd_start,
                    'tics_start':   tics_start,
                    'tjd_stop':     tjd_stop,
                    'tics_stop':    tics_stop,
                    'tstart':       tstart,
                    'tstop':        tstop,
                    'gti':          gti,
                    'quality':      quality}

        # Close FITS file
        fits.close()

        # Return information
        return oad_dict

    def _get_hkd_list(self):
        """
        Get HKD list

        Returns
        -------
        hkd_list : list of dict
            HKD database
        """
        # Initialse HKD list
        hkd_list = []

        # Get available phases
        phases = self._get_phases()

        # Loop over phases
        for phase in phases:

            # Get available VPs
            vps = glob.glob('%s/vp*' % phase)
            vps.sort()

            # Loop over VPs
            for vp in vps:

                # Get HKDs for VP
                hkds = self._get_hkds_for_vp(vp)

                # Append HKDs to list
                hkd_list.extend(hkds)

        # Sort HK list by start day
        hkd_list = self._sorted(hkd_list)

        # Return HKD list
        return hkd_list

    def _get_hkds_for_vp(self, vp):
        """
        Build HKD database entries for Viewing Period

        Parameters
        ----------
        vp : str
            VP path

        Returns
        -------
        hkd_list : list of dict
            HKD database entries
        """
        # Get archive name and expanded path
        archive      = self['archive'].string()
        archive_path = gammalib.expand_env(archive)

        # Get VP name
        vpname = os.path.basename(vp)

        # Initialse HKD list
        hkd_list = []

        # Get HKD files
        hkds = glob.glob('%s/*hkd.fits*' % vp)

        # Loop over all HKD files
        for hkd in hkds:

            # Set HKD filename
            hkdname = archive+'/'+self._get_relpath(hkd, archive_path).strip('.gz')

            # Try opening file
            try:
                fits = gammalib.GFits(hkdname)
            except:
                self._log_string(gammalib.NORMAL, '*** Unable to open HKD file "%s". '
                                 'Skip file.' % (hkdname))
                continue

            # Get information from HKD file
            hkd_dict = self._get_hkd_info(hkdname)

            # Skip file if it has a negative quality
            if hkd_dict['quality'] < 0:
                self._log_string(gammalib.NORMAL, '*** HKD file "%s" has quality %d. '
                                 'Skip file.' % (hkdname, hkd_dict['quality']))
                continue

            # Set VP
            hkd_dict['vp'] = vpname

            # Append dictionary
            hkd_list.append(hkd_dict)

        # Sort HKD list by start day
        hkd_list = self._sorted(hkd_list)

        # Return HKD list
        return hkd_list

    def _get_hkd_info(self, hkdname):
        """
        Build HKD database entry from HKD file

        Parameters
        ----------
        hkdname : str
            HKD filename

        Returns
        -------
        hkd_dict : dict
            HKD database entry
        """
        # Open OAD file
        fits = gammalib.GFits(hkdname)

        # Get superpacket table
        records  = fits.table(1)
        nrecords = records.nrows()

        # Extract attributes for HKD file
        obs_id   = records.real('OBS_ID')
        object   = records.string('OBJECT')
        if records.has_card('GLON_SCZ'):
            glon_scz = records.real('GLON_SCZ')
        else:
            glon_scz = -1000.0
        if records.has_card('GLAT_SCZ'):
            glat_scz = records.real('GLAT_SCZ')
        else:
            glat_scz = -1000.0
        quality = records.integer('DSD_QUA')
        dsd_id  = '%s-%s-%d' % (records.string('DSD_SB'), records.string('DSTID'),
                                records.integer('DSD_NO'))

        # Get time columns
        tjd  = records['TJD']
        tics = records['TICS']

        # Get start and end time. Catch any errors since files may be corrupt.
        # In case of errors, use times from file header.
        try:
            tjd_start  = tjd[0]
            tics_start = tics[0]
            tjd_stop   = tjd[nrecords-1]
            tics_stop  = tics[nrecords-1]
        except:
            self._log_string(gammalib.NORMAL, '*** Corrupt HKD file "%s" encountered. '
                             'Set HKD quality to -200.' % (hkdname))
            tjd_start  = records.integer('VISDAY')
            tics_start = records.integer('VISTIM')
            tjd_stop   = records.integer('VIEDAY')
            tics_stop  = records.integer('VIETIM')
            quality    = -200

        # Compute GTime of first and last entry
        tstart = gammalib.com_time(tjd_start, tics_start)
        tstop  = gammalib.com_time(tjd_stop,  tics_stop)

        # Set GTI for HKD
        gti = gammalib.GGti(tstart, tstop)

        # Build dictionary
        hkd_dict = {'filename':   hkdname,
                    'records':    nrecords,
                    'obs_id':     obs_id,
                    'dsd_id':     dsd_id,
                    'object':     object,
                    'glon_scz':   glon_scz,
                    'glat_scz':   glat_scz,
                    'date_start': tstart.utc(),
                    'date_stop':  tstop.utc(),
                    'tjd_start':  tjd_start,
                    'tics_start': tics_start,
                    'tjd_stop':   tjd_stop,
                    'tics_stop':  tics_stop,
                    'tstart':     tstart,
                    'tstop':      tstop,
                    'gti':        gti,
                    'quality':    quality}

        # Close FITS file
        fits.close()

        # Return information
        return hkd_dict

    def _get_tim_list(self):
        """
        Get TIM list

        Returns
        -------
        tim_list : list of dict
            TIM database
        """
        # Initialse TIM list
        tim_list = []

        # Get available phases
        phases = self._get_phases()

        # Loop over phases
        for phase in phases:

            # Get available VPs
            vps = glob.glob('%s/vp*' % phase)
            vps.sort()

            # Loop over VPs
            for vp in vps:

                # Get TIMs for VP
                tims = self._get_tims_for_vp(vp)

                # Append TIMs to list
                tim_list.extend(tims)

        # Sort TIM list by start day
        tim_list = self._sorted(tim_list)

        # Return TIM list
        return tim_list

    def _get_tims_for_vp(self, vp):
        """
        Build TIM database entries for Viewing Period

        Parameters
        ----------
        vp : str
            VP path

        Returns
        -------
        tim_list : list of dict
            TIM database entries
        """
        # Get archive name and expanded path
        archive      = self['archive'].string()
        archive_path = gammalib.expand_env(archive)

        # Get VP name
        vpname = os.path.basename(vp)

        # Initialse TIM list
        tim_list = []

        # Get TIM files
        tims = glob.glob('%s/*tim.fits*' % vp)

        # Loop over all TIM files
        for tim in tims:

            # Set TIM filename
            timname = archive+'/'+self._get_relpath(tim, archive_path).strip('.gz')

            # Try opening file
            try:
                fits = gammalib.GFits(timname)
            except:
                self._log_string(gammalib.NORMAL, '*** Unable to open TIM file "%s". '
                                 'Skip file.' % (timname))
                continue

            # Get information from TIM file
            tim_dict = self._get_tim_info(timname)

            # Skip file if it has a negative quality
            if tim_dict['quality'] < 0:
                self._log_string(gammalib.NORMAL, '*** TIM file "%s" has quality %d. '
                                 'Skip file.' % (timname, tim_dict['quality']))
                continue

            # Set VP
            tim_dict['vp'] = vpname

            # Append dictionary
            tim_list.append(tim_dict)

        # Sort TIM list by start day
        tim_list = self._sorted(tim_list)

        # Return TIM list
        return tim_list

    def _get_tim_info(self, timname):
        """
        Build TIM database entry from TIM file

        Parameters
        ----------
        timname : str
            TIM filename

        Returns
        -------
        tim_dict : dict
            TIM database entry
        """
        # Open TIM file
        fits = gammalib.GFits(timname)

        # Get superpacket table
        records  = fits.table(1)
        nrecords = records.nrows()

        # Extract attributes for TIM file
        obs_id   = records.real('OBS_ID')
        object   = records.string('OBJECT')
        visday   = records.integer('VISDAY')
        vistim   = records.integer('VISTIM')
        vieday   = records.integer('VIEDAY')
        vietim   = records.integer('VIETIM')
        glon_scz = records.real('GLON_SCZ')
        glat_scz = records.real('GLAT_SCZ')
        quality  = records.integer('DSD_QUA')
        dsd_id   = '%s-%s-%d' % (records.string('DSD_SB'), records.string('DSTID'),
                                 records.integer('DSD_NO'))

        # Get time columns
        start_tjd = records['START_TJD']
        start_tic = records['START_TIC']
        end_tjd   = records['END_TJD']
        end_tic   = records['END_TIC']

        # Get start and end time. Catch any errors since files may be corrupt.
        # In case of errors, use times from file header.
        try:
            tjd_start  = start_tjd[0]
            tics_start = start_tic[0]
            tjd_stop   = end_tjd[nrecords-1]
            tics_stop  = end_tic[nrecords-1]
        except:
            self._log_string(gammalib.NORMAL, '*** Corrupt TIM file "%s" encountered. '
                             'Set TIM quality to -200.' % (evpname))
            tjd_start  = visday
            tics_start = vistim
            tjd_stop   = vieday
            tics_stop  = vietim
            quality    = -200

        # Close FITS file
        fits.close()

        # Compute GTime of first and last superpacket
        tstart = gammalib.com_time(tjd_start, tics_start)
        tstop  = gammalib.com_time(tjd_stop, tics_stop)

        # Set GTI for TIM
        gti = gammalib.GGti(tstart, tstop)

        # Get telapse and ontime
        tim     = gammalib.GCOMTim(timname)
        telapse = tim.gti().telapse()
        ontime  = tim.gti().ontime()

        # Build dictionary
        tim_dict = {'filename':   timname,
                    'intervals':  nrecords,
                    'obs_id':     obs_id,
                    'dsd_id':     dsd_id,
                    'object':     object,
                    'glon_scz':   glon_scz,
                    'glat_scz':   glat_scz,
                    'date_start': tstart.utc(),
                    'date_stop':  tstop.utc(),
                    'tjd_start':  tjd_start,
                    'tics_start': tics_start,
                    'tjd_stop':   tjd_stop,
                    'tics_stop':  tics_stop,
                    'visday':     visday,
                    'vistim':     vistim,
                    'vieday':     vieday,
                    'vietim':     vietim,
                    'tstart':     tstart,
                    'tstop':      tstop,
                    'telapse':    telapse,
                    'ontime':     ontime,
                    'gti':        gti,
                    'quality':    quality}

        # Return information
        return tim_dict

    def _get_xml_list(self, evp_list, tim_list, oad_list, hkd_list):
        """
        Build XML database from EVP, TIM, OAD and HKD database

        Parameters
        ----------
        evp_list : list of dict
            EVP database
        tim_list : list of dict
            TIM database
        oad_list : list of dict
            OAD database
        hkd_list : list of dict
            HKD database

        Returns
        -------
        xml_list : list of dict
            XML database
        """
        # Initialse XML list
        xml_list = []

        # Get filtered TIM list
        filtered_tim_list = self._filter_tim_list(tim_list)

        # Get available phases
        phases = self._get_phases()

        # Loop over phases
        for phase in phases:

            # Get available VPs
            vps = glob.glob('%s/vp*' % phase)
            vps.sort()

            # Loop over VPs
            for vp in vps:

                # Get XML for VP
                xml = self._get_xml_for_vp(vp, evp_list, filtered_tim_list, oad_list, hkd_list)

                # Append XML to list
                if xml is not None:
                    xml_list.append(xml)

        # Sort XML list by start day
        xml_list = self._sorted(xml_list)

        # Return XML list
        return xml_list

    def _filter_tim_list(self, tim_list):
        """
        Removes all invalid TIMs from TIM database

        For a given VP only the TIMs with the maximum quality will be kept. This
        is needed since some of the TIM dataset cover periods well beyond a given
        VP, and these TIMs are marked with data quality 135. On the other hand,
        there are VPs with valid TIMs of data quality 135, hence we cannot simply
        filter on data quality.

        Parameters
        ----------
        tim_list : list of dict
            TIM database

        Returns
        -------
        filtered_tim_list : list of dict
            Filtered TIM database
        """
        # Initialise filtered TIM list
        filtered_tim_list = []

        # Get available phases
        phases = self._get_phases()

        # Loop over phases
        for phase in phases:

            # Get available VPs
            vps = glob.glob('%s/vp*' % phase)
            vps.sort()

            # Loop over VPs
            for vp in vps:

                # Get VP name
                vpname = os.path.basename(vp)

                # Get all TIMs for VP
                tims = list(filter(lambda k: k['vp'] == vpname, tim_list))

                # Continue if list is empty
                if len(tims) < 1:
                    continue

                # Get maximum quality of all TIMs in VP
                quality = max([tim['quality'] for tim in tims])

                # Get all TIMs for VP with maximum quality
                best_tims = list(filter(lambda k: k['quality'] == quality, tims))

                # Append these TIMs to the filtered list
                filtered_tim_list.extend(best_tims)

        # Return
        return filtered_tim_list

    def _get_xml_for_vp(self, vp, evp_list, tim_list, oad_list, hkd_list):
        """
        Build XML database entry for Viewing Period

        Parameters
        ----------
        vp : str
            Viewing Period directory
        evp_list : list of dict
            EVP database
        tim_list : list of dict
            TIM database
        oad_list : list of dict
            OAD database
        hkd_list : list of dict
            HKD database

        Returns
        -------
        xml_dict : dictionary
            XML database entry
        """
        # Get database directory
        path = self['dbase'].string()

        # Set VP name
        vpname  = os.path.basename(vp)

        # Get EVP for VP
        evp = self._get_evp(vp, evp_list)
        if evp is None:
            if vpname == 'vp0225_0' or \
               vpname == 'vp0303_5' or \
               vpname == 'vp0308_3' or \
               vpname == 'vp0311_3' or \
               vpname == 'vp0617_2' or \
               vpname == 'vp0617_4' or \
               vpname == 'vp0617_6' or \
               vpname == 'vp0619_5':
                self._log_string(gammalib.NORMAL, '--- No EVP file found for "%s due '
                                 'to Reboost". Skip VP.' % (vp))
            elif vpname == 'vp0229_3':
                self._log_string(gammalib.NORMAL, '--- No EVP file found for "%s due '
                                 'to Perseids". Skip VP.' % (vp))
            elif vpname == 'vp0833_7':
                self._log_string(gammalib.NORMAL, '--- No EVP file found for "%s due '
                                 'to Leonid Shower". Skip VP.' % (vp))
            else:
                self._log_string(gammalib.NORMAL, '*** No EVP file found for "%s". '
                                 'Skip VP.' % (vp))
            return None

        # Get TIM for EVP
        tim = self._get_tim(evp, tim_list)
        if tim is None:
            self._log_string(gammalib.NORMAL, '*** No TIM file found for "%s". '
                             'Skip VP.' % (vp))
            return None

        # Get OADs for EVP
        oads = self._get_oads(evp, oad_list)
        if len(oads) == 0:
            self._log_string(gammalib.NORMAL, '*** No OAD files found for "%s". '
                             'Skip VP.' % (vp))
            return None

        # Get HKDs for EVP
        hkds = self._get_hkds(evp, hkd_list)
        if len(hkds) == 0:
            self._log_string(gammalib.NORMAL, '*** No HKD files found for "%s". '
                             'Skip VP.' % (vp))
            return None

        # Remove OADs that are not covered by HKDs
        oads_selected = []
        for oad in oads:
            tjd_start = oad['tjd_start']
            tjd_stop  = oad['tjd_stop']
            match     = False
            for hkd in hkds:
                if tjd_start == hkd['tjd_start'] and \
                   tjd_stop  == hkd['tjd_stop']:
                   match = True
                   break
            if not match:
                self._log_string(gammalib.NORMAL, '*** No HKD file found for TJD '
                                 '%d for "%s". Remove corresponding OAD.' %
                                 (tjd_start, vp))
            else:
                oads_selected.append(oad)
        oads = oads_selected

        # Get oadnames and oad qualities
        oadnames     = [oad['filename'] for oad in oads]
        oadqualities = [oad['quality']  for oad in oads]

        # Get hkdnames and hkd qualities
        hkdnames     = [hkd['filename'] for hkd in hkds]
        hkdqualities = [hkd['quality']  for hkd in hkds]

        # Build dictionary
        xml_dict = {'vp':          vpname,
                    'evpname':     evp['filename'],
                    'timname':     tim['filename'],
                    'oadnames':    oadnames,
                    'hkdnames':    hkdnames,
                    'xmlname':     '%s/xml/%s.xml' % (path, vpname),
                    'nevents':     evp['nevents'],
                    'obs_id':      evp['obs_id'],
                    'object':      evp['object'],
                    'glon_scz':    evp['glon_scz'],
                    'glat_scz':    evp['glat_scz'],
                    'date_start':  tim['date_start'],
                    'date_stop':   tim['date_stop'],
                    'tjd_start':   tim['tjd_start'],
                    'tics_start':  tim['tics_start'],
                    'tjd_stop':    tim['tjd_stop'],
                    'tics_stop':   tim['tics_stop'],
                    'evp_tstart':  evp['tstart'],
                    'evp_tstop':   evp['tstop'],
                    'visday':      tim['visday'],
                    'vistim':      tim['vistim'],
                    'vieday':      tim['vieday'],
                    'vietim':      tim['vietim'],
                    'telapse':     tim['telapse'],
                    'ontime':      tim['ontime'],
                    'evp_verno':   evp['verno'],
                    'evp_quality': evp['quality'],
                    'tim_quality': tim['quality'],
                    'oad_quality': oadqualities,
                    'hkd_quality': hkdqualities}

        # Return XML list
        return xml_dict

    def _get_evp(self, vp, evp_list):
        """
        Get EVP database entry for Viewing Period

        Parameters
        ----------
        vp : str
            Viewing Period directory
        evp_list : list of dict
            EVP database

        Returns
        -------
        best_evp : dict
            EVP database entry
        """
        # Get quality flag
        qf = self['quality'].integer()

        # Get VP name
        vpname = os.path.basename(vp)

        # Get all EVPs for VP
        evps = list(filter(lambda k: k['vp'] == vpname, evp_list))

        # Initialise best EVP
        best_verno_evp   = None
        best_quality_evp = None
        best_verno       = 0
        best_quality     = 0

        # Initialise choosen flag
        choosen = False

        # First find EVP with target QF and highest version number
        for evp in evps:
            if evp['quality'] == qf:
                if evp['verno'] > best_verno:
                    choosen    = True
                    best_evp   = evp
                    best_verno = evp['verno']

        # Continue only if no EVP was found
        if not choosen:

            # Find EVP with highest version number
            for evp in evps:
                if evp['verno'] > best_verno:
                    best_verno     = evp['verno']
                    best_verno_evp = evp

            # Find EVP with best quality
            for evp in evps:
                if evp['quality'] > best_quality:
                    best_quality     = evp['quality']
                    best_quality_evp = evp

            # Initialise best EVP with best quality EVP
            best_evp = best_quality_evp

            # If highest version number is not the EVP with the best quality then
            # select the EVP with the highest version number
            if best_verno_evp != None and best_quality_evp != None:
                if best_verno_evp['filename'] != best_quality_evp['filename']:
                    if best_quality_evp['verno'] < best_verno:

                        # Signal issue
                        self._log_string(gammalib.NORMAL, '*** EVP "%s" with best '
                                         'quality has version %d which is not the '
                                         'highest version %d.' %
                                         (best_quality_evp['filename'],
                                          best_quality_evp['verno'],
                                          best_verno))

                        # If there are several EVPs with highest version then select
                        # the one with the best quality
                        evps = list(filter(lambda k: k['vp'] == vpname and k['verno'] == best_verno, evp_list))
                        if len(evps) > 1:
                            self._log_string(gammalib.NORMAL, '--- Several EVPs for %s '
                                             'have highest version' % (vpname))
                            best = 0
                            for evp in evps:
                                if evp['quality'] > best:
                                    best     = evp['quality']
                                    best_evp = evp
                            for evp in evps:
                                if best_evp['filename'] == evp['filename']:
                                    choose  = '+'
                                    choosen = True
                                else:
                                    choose = ' '
                                self._log_string(gammalib.NORMAL, '  %s "%s" (verno=%d, '
                                                 'quality=%d)' %
                                                 (choose,
                                                  evp['filename'],
                                                  evp['verno'],
                                                  evp['quality']))

                        # ... otherwise use the EVP with the best version number
                        else:
                            best_evp = best_verno_evp
                            choosen  = True
                            self._log_string(gammalib.NORMAL, '--- Choose "%s" (verno=%d, '
                                             'quality=%d)' %
                                             (best_evp['filename'],
                                              best_evp['verno'],
                                              best_evp['quality']))

            # Continue only if no file was choosen
            if not choosen:

                # If there are several EVPs with the best quality then select the one
                # with the highest version
                evps = list(filter(lambda k: k['vp'] == vpname and k['quality'] == best_quality, evp_list))
                if len(evps) > 1:
                    self._log_string(gammalib.NORMAL, '--- Several EVPs for %s have '
                                     'best quality' % (vpname))
                    best = 0
                    for evp in evps:
                        if evp['verno'] > best:
                            best     = evp['verno']
                            best_evp = evp
                    for evp in evps:
                        if best_evp['filename'] == evp['filename']:
                            choose = '+'
                        else:
                            choose = ' '
                        self._log_string(gammalib.NORMAL, '  %s "%s" (verno=%d, '
                                         'quality=%d)' %
                                         (choose,
                                          evp['filename'],
                                          evp['verno'],
                                          evp['quality']))

        # Return best EVP
        return best_evp

    def _get_tim(self, evp, tim_list):
        """
        Get TIM database entry for EVP database entry

        Parameters
        ----------
        evp : dict
            EVP dictionary
        tim_list : list of dict
            TIM database

        Returns
        -------
        best_tim : dict
            TIM database entry for EVP database entry
        """
        # Initialise maximum overlap
        max_overlap = 0.0
        min_dtabs   = 1.0e30
        best_tim    = None

        # Loop over all TIMs
        for tim in tim_list:

            # Compute time difference to EVP start time and stop time
            dtstart = evp['tstart'] - tim['tstart']
            dtstop  = evp['tstop']  - tim['tstop']

            # Compute absolute time difference
            dtabs = abs(dtstart) + abs(dtstop)

            # Keep smallest time difference
            if dtabs < min_dtabs:
                min_dtabs = dtabs
                best_tim  = tim

        # Return best TIM
        return best_tim

    def _get_oads(self, evp, oad_list):
        """
        Get OAD database entry for EVP database entry

        Parameters
        ----------
        evp : dict
            EVP dictionary
        oad_list : list of dict
            OAD database

        Returns
        -------
        oads : list of dict
            OAD database entries for EVP database entry
        """
        # Initialise list of OADs
        oads = []

        # Loop over all OADs
        for oad in oad_list:

            # Compute overlap between EVP and OAD
            overlap = evp['gti'].overlap(oad['tstart'], oad['tstop'])

            # If we have overlap then check if the file already exists
            if overlap > 0.0:

                # First check if OAD is already in OAD list. If this is the case,
                # and the current OAD corresponds to the VP of the EVP, then replace
                # the existing OAD by the current OAD since we prefer having the
                # OAD in the same folder as the EVP. Otherwise we keep the exiting
                # OAD
                exists = False
                for i, oad_exist in enumerate(oads):
                    if oad_exist['dsd_id'] == oad['dsd_id']:
                        exists = True
                        if oad['vp'] == evp['vp']:
                            oads[i] = oad
                            break
                        else:
                            break

                # If the OAD does not yet exist then append it now
                if not exists:
                    oads.append(oad)

        # Sort OAD list by start day
        oads = self._sorted(oads)

        # Return OADs
        return oads

    def _get_hkds(self, evp, hkd_list):
        """
        Get HKD database entry for EVP database entry

        Parameters
        ----------
        evp : dict
            EVP dictionary
        hkd_list : list of dict
            HKD database

        Returns
        -------
        hkds : list of dict
            HKD database entries for EVP database entry
        """
        # Initialise list of HKDs
        hkds = []

        # Loop over all HKDs
        for hkd in hkd_list:

            # Compute overlap between EVP and HKD
            overlap = evp['gti'].overlap(hkd['tstart'], hkd['tstop'])

            # If we have overlap then check if the file already exists
            if overlap > 0.0:

                # First check if HKD is already in HKD list. If this is the case,
                # and the current HKD corresponds to the VP of the EVP, then replace
                # the existing HKD by the current HKD since we prefer having the
                # HKD in the same folder as the EVP. Otherwise we keep the existing
                # HKD
                exists = False
                for i, hkd_exist in enumerate(hkds):
                    if hkd_exist['dsd_id'] == hkd['dsd_id']:
                        exists = True
                        if hkd['vp'] == evp['vp']:
                            hkds[i] = hkd
                            break
                        else:
                            break

                # If the HKD does not yet exist then append it now
                if not exists:
                    hkds.append(hkd)

        # Sort HKD list by start day
        hkds = self._sorted(hkds)

        # Return HKDs
        return hkds

    def _check_xml_list(self, xml_list):
        """
        Check XML list against COMPASS database

        Returns
        -------
        xml_list : list of dict
            XML database with mode and status information
        """
        # Load COMPASS database
        database = self._get_compass_database()

        # Loop over XML list
        for xml_dict in xml_list:

            # Add mode and status to XML list
            xml_dict['mode']   = '???'
            xml_dict['status'] = 'unknown'

            # Continue only if there are entries
            if len(database) > 0:

                # Find VP in COMPASS database
                compass = None
                for entry in database:
                    if entry['vp'] == xml_dict['vp']:
                        compass = entry
                        break

                # Signal if VP was not found in COMPASS database
                if compass == None:
                    xml_dict['status'] = 'not found'
                    self._log_string(gammalib.NORMAL, '--- VP "%s" not found in COMPASS '
                                     'database.' % (xml_dict['vp']))
                    continue

                # Set mode and status
                xml_dict['mode']   = compass['mode']
                xml_dict['status'] = 'found'

                # Check EVP filename
                evp = os.path.basename(xml_dict['evpname'])[1:].rstrip('_evp.fits')
                if evp != compass['evp']:
                    self._log_string(gammalib.NORMAL, '--- VP "%s": EVP "%s" in COMPASS '
                                     'database differs from EVP "%s" in XML database.' %
                                     (xml_dict['vp'], compass['evp'], evp))

                # Check TIM filename
                tim = os.path.basename(xml_dict['timname'])[1:].rstrip('_tim.fits')
                if tim != compass['tim']:
                    self._log_string(gammalib.NORMAL, '--- VP "%s": TIM "%s" in COMPASS '
                                     'database differs from TIM "%s" in XML database.' %
                                     (xml_dict['vp'], compass['tim'], tim))

        # Return
        return xml_list

    def _get_compass_database(self):
        """
        Get COMPASS database

        Returns
        -------
        database : list of dict
            COMPASS database
        """
        # Initialise database
        database = []

        # Continue only if database was specified
        if self['refdata'].is_valid():

            # Read ASCII file
            f = open(gammalib.expand_env(self['refdata'].filename().url()), 'r')
            for line in f:

                # Skip comment lines
                if line[0:1] == 'C':
                    continue

                # Format VP
                vp = ('vp%06.1f' % (float(line[1:6]))).replace('.', '_')

                # Extract fields
                record = {'vp':   vp,
                          'evp':  gammalib.strip_whitespace(line[13:19]),
                          'tim':  gammalib.strip_whitespace(line[21:27]),
                          'oad':  gammalib.strip_whitespace(line[29:35]),
                          'mode': gammalib.strip_whitespace(line[66:70])}

                # Append record to database
                database.append(record)

        # Return database
        return database

    def _create_dbase_dirs(self):
        """
        Create database directory structure
        """
        # Get database directory
        path = self['dbase'].string()

        # Create database directory
        try:
            os.makedirs(gammalib.expand_env(path))
        except OSError:
            self._log_string(gammalib.NORMAL, '--- Directory "%s" exists' % path)
        try:
            os.makedirs('%s/xml' % gammalib.expand_env(path))
        except OSError:
            self._log_string(gammalib.NORMAL, '--- Directory "%s/xml" exists' % path)
        try:
            os.makedirs('%s/tim' % gammalib.expand_env(path))
        except OSError:
            self._log_string(gammalib.NORMAL, '--- Directory "%s/tim" exists' % path)

        # Return
        return

    def _build_tim_file(self, xml_dict):
        """
        Build TIM file for XML database entry

        Parameters
        ----------
        xml_dict : dict
            XML database entry
        """
        # Get database directory
        path = self['dbase'].string()

        # Get TIM from database
        tim = gammalib.GCOMTim(xml_dict['timname'])

        # Extract GTIs
        gti_heasarc = tim.gti()

        # Get tstart and tstop from first and last event
        tstart = xml_dict['evp_tstart']
        tstop  = xml_dict['evp_tstop']

        # Kluge for VP vp0411_1 which overlaps with vp0411_5
        if xml_dict['vp'] == 'vp0411_1':
            tstop = gammalib.com_time(9769, 444782835)
            self._log_string(gammalib.NORMAL, '*** Stop date for "%s" changed from %s '
                             'to %s.' %
                             (xml_dict['vp'], xml_dict['evp_tstop'].utc(), tstop.utc()))

        # Initialise GTI
        gti = gammalib.GGti()

        # Select GTIs
        search_start = True
        for i in range(len(gti_heasarc)):

            # Case A: search for start
            if search_start:

                # If event start time is after GTI stop time then skip this GTI
                if tstart > gti_heasarc.tstop(i):
                    continue

                # ... otherwise stop searching for start and set first GTI
                else:
                    # Stop searching for start
                    search_start = False

                    # Get first GTI
                    if tstart < gti_heasarc.tstart(i):
                        t_start = gti_heasarc.tstart(i)
                        t_stop  = gti_heasarc.tstop(i)
                    else:
                        t_start = tstart
                        t_stop  = gti_heasarc.tstop(i)

                    # Append GTI
                    gti.append(t_start, t_stop)

            # Case B: search for stop
            else:

                # If event stop time after GTI stop time then and this GTI
                if tstop > gti_heasarc.tstop(i):
                    t_start = gti_heasarc.tstart(i)
                    t_stop  = gti_heasarc.tstop(i)
                    gti.append(t_start, t_stop)

                # ... otherwise, if the event stop time is before the GTI start
                # time then stop appending GTIs
                elif tstop < gti_heasarc.tstart(i):
                    break

                # ... otherwise append last GTI
                else:
                    t_start = gti_heasarc.tstart(i)
                    t_stop  = tstop
                    gti.append(t_start, t_stop)
                    break

        # Set new GTIs for TIM
        tim.gti(gti)

        # Set new TIM filename
        xml_dict['timname'] = '%s/tim/%s_tim.fits' % (path, xml_dict['vp'])

        # Save TIM
        tim.save(xml_dict['timname'], True)

        # Return
        return

    def _build_xml_file(self, xml_dict):
        """
        Build observation definition XML file for XML database entry

        Parameters
        ----------
        xml_dict : dict
            XML database entry
        """
        # Create XML object
        xml = gammalib.GXml()

        # Append observation list
        obslist = xml.append('observation_list title="observation list"')

        # Append observation
        obs = obslist.append('observation name="%s" id="%s" instrument="COM"' % \
                             (xml_dict['object'], xml_dict['vp']))

        # Append EVP file
        obs.append('parameter name="EVP" file="%s"' % (xml_dict['evpname']))

        # Append TIM file
        obs.append('parameter name="TIM" file="%s"' % (xml_dict['timname']))

        # Append OAD files
        for oadname in xml_dict['oadnames']:
            obs.append('parameter name="OAD" file="%s"' % (oadname))

        # Append HKD files
        for hkdname in xml_dict['hkdnames']:
            obs.append('parameter name="HKD" file="%s"' % (hkdname))

        # Save XML file
        xml.save(xml_dict['xmlname'])

        # Return
        return

    def _check_database(self, xml_list):
        """
        Check database

        Parameters
        ----------
        xml_list : list of dict
            XML database
        """
        # Determine number of database entries
        nrows = len(xml_list)

        # Check that stop dates are after start dates
        for i in range(nrows):
            tstart = gammalib.GTime()
            tstop  = gammalib.GTime()
            tstart.utc(xml_list[i]['date_start'])
            tstop.utc(xml_list[i]['date_stop'])
            if tstart > tstop:
                self._log_string(gammalib.NORMAL, '*** Start date %s of VP %s later '
                                 'than stop date %s.' %
                                 (xml_list[i]['date_start'],
                                  xml_list[i]['vp'],
                                  xml_list[i]['date_stop']))

        # Check that stop dates of VP is before start date of next VP. If this is
        # not the case, events may be double counted.
        for i in range(nrows-1):
            tstart = gammalib.GTime()
            tstop  = gammalib.GTime()
            tstop.utc(xml_list[i]['date_stop'])
            tstart.utc(xml_list[i+1]['date_start'])
            if tstop >= tstart:
                self._log_string(gammalib.NORMAL, '*** Stop date %s of VP %s later '
                                 'than start date %s of VP %s.' %
                                 (xml_list[i]['date_stop'],
                                  xml_list[i]['vp'],
                                  xml_list[i+1]['date_start'],
                                  xml_list[i+1]['vp']))

        # Return
        return

    def _build_database_file(self, evp_list, oad_list, hkd_list, tim_list, xml_list):
        """
        Build database file

        Parameters
        ----------
        evp_list : list of dict
            EVP database
        oad_list : list of dict
            OAD database
        hkd_list : list of dict
            HKD database
        tim_list : list of dict
            TIM database
        xml_list : list of dict
            XML database
        """
        # Setup filename
        filename = '%s/dbase.fits' % (self['dbase'].string())

        # Create FITS file
        fits = gammalib.GFits()

        # Get tables
        evp_table = self._build_evp_table(evp_list)
        oad_table = self._build_oad_table(oad_list)
        hkd_table = self._build_hkd_table(hkd_list)
        tim_table = self._build_tim_table(tim_list)
        xml_table = self._build_xml_table(xml_list)

        # Append tables
        fits.append(evp_table)
        fits.append(oad_table)
        fits.append(hkd_table)
        fits.append(tim_table)
        fits.append(xml_table)

        # Save database
        fits.saveto(filename, True)

        # Close FITS file
        fits.close()

        # Return
        return

    def _build_evp_table(self, evp_list):
        """
        Build EVP FITS table

        Parameters
        ----------
        evp_list : list of dict
            EVP database

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            EVP FITS table
        """
        # Determine number of rows in table
        nrows = len(evp_list)

        # Create FITS table columns
        col_vp         = gammalib.GFitsTableStringCol('VP',         nrows, 9)
        col_object     = gammalib.GFitsTableStringCol('OBJECT',     nrows, 20)
        col_obs_id     = gammalib.GFitsTableStringCol('OBS_ID',     nrows, 9)
        col_dsd_id     = gammalib.GFitsTableStringCol('DSD_ID',     nrows, 15)
        col_glon_scz   = gammalib.GFitsTableFloatCol('GLON_SCZ',    nrows)
        col_glat_scz   = gammalib.GFitsTableFloatCol('GLAT_SCZ',    nrows)
        col_date_start = gammalib.GFitsTableStringCol('DATE_START', nrows, 20)
        col_date_stop  = gammalib.GFitsTableStringCol('DATE_STOP',  nrows, 20)
        col_tjd_start  = gammalib.GFitsTableShortCol('TJD_START',   nrows)
        col_tics_start = gammalib.GFitsTableLongCol('TICS_START',   nrows)
        col_tjd_stop   = gammalib.GFitsTableShortCol('TJD_STOP',    nrows)
        col_tics_stop  = gammalib.GFitsTableLongCol('TICS_STOP',    nrows)
        col_visday     = gammalib.GFitsTableShortCol('VISDAY',      nrows)
        col_vistim     = gammalib.GFitsTableLongCol('VISTIM',       nrows)
        col_vieday     = gammalib.GFitsTableShortCol('VIEDAY',      nrows)
        col_vietim     = gammalib.GFitsTableLongCol('VIETIM',       nrows)
        col_verno      = gammalib.GFitsTableShortCol('VERNO',       nrows)
        col_quality    = gammalib.GFitsTableShortCol('QUALITY',     nrows)
        col_evp        = gammalib.GFitsTableStringCol('EVP',        nrows, 60)
        col_nevents    = gammalib.GFitsTableLongCol('NEVENTS',      nrows)

        # Set column units
        col_glon_scz.unit('deg')
        col_glat_scz.unit('deg')
        col_tjd_start.unit('days')
        col_tics_start.unit('0.000125 s')
        col_tjd_stop.unit('days')
        col_tics_stop.unit('0.000125 s')
        col_visday.unit('days')
        col_vistim.unit('0.000125 s')
        col_vieday.unit('days')
        col_vietim.unit('0.000125 s')

        # Fill columns
        for i, evp_record in enumerate(evp_list):
            col_vp[i]         = evp_record['vp']
            col_object[i]     = evp_record['object']
            col_obs_id[i]     = '%.1f' % evp_record['obs_id']
            col_dsd_id[i]     = evp_record['dsd_id']
            col_glon_scz[i]   = evp_record['glon_scz']
            col_glat_scz[i]   = evp_record['glat_scz']
            col_date_start[i] = evp_record['date_start']
            col_date_stop[i]  = evp_record['date_stop']
            col_tjd_start[i]  = evp_record['tjd_start']
            col_tics_start[i] = evp_record['tics_start']
            col_tjd_stop[i]   = evp_record['tjd_stop']
            col_tics_stop[i]  = evp_record['tics_stop']
            col_visday[i]     = evp_record['visday']
            col_vistim[i]     = evp_record['vistim']
            col_vieday[i]     = evp_record['vieday']
            col_vietim[i]     = evp_record['vietim']
            col_verno[i]      = evp_record['verno']
            col_quality[i]    = evp_record['quality']
            col_evp[i]        = evp_record['filename']
            col_nevents[i]    = evp_record['nevents']

        # Create FITS table
        table = gammalib.GFitsBinTable()

        # Append columns to table
        table.append(col_vp)
        table.append(col_object)
        table.append(col_obs_id)
        table.append(col_dsd_id)
        table.append(col_glon_scz)
        table.append(col_glat_scz)
        table.append(col_date_start)
        table.append(col_date_stop)
        table.append(col_tjd_start)
        table.append(col_tics_start)
        table.append(col_tjd_stop)
        table.append(col_tics_stop)
        table.append(col_visday)
        table.append(col_vistim)
        table.append(col_vieday)
        table.append(col_vietim)
        table.append(col_verno)
        table.append(col_quality)
        table.append(col_evp)
        table.append(col_nevents)

        # Set extension name
        table.extname('EVP')

        # Return table
        return table

    def _build_oad_table(self, oad_list):
        """
        Build OAD FITS table

        Parameters
        ----------
        oad_list : list of dict
            OAD database

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            OAD FITS table
        """
        # Determine number of rows in table
        nrows = len(oad_list)

        # Create FITS table columns
        col_vp           = gammalib.GFitsTableStringCol('VP',         nrows, 9)
        col_object       = gammalib.GFitsTableStringCol('OBJECT',     nrows, 20)
        col_obs_id       = gammalib.GFitsTableStringCol('OBS_ID',     nrows, 9)
        col_dsd_id       = gammalib.GFitsTableStringCol('DSD_ID',     nrows, 15)
        col_glon_scz     = gammalib.GFitsTableFloatCol('GLON_SCZ',    nrows)
        col_glat_scz     = gammalib.GFitsTableFloatCol('GLAT_SCZ',    nrows)
        col_date_start   = gammalib.GFitsTableStringCol('DATE_START', nrows, 20)
        col_date_stop    = gammalib.GFitsTableStringCol('DATE_STOP',  nrows, 20)
        col_tjd_start    = gammalib.GFitsTableShortCol('TJD_START',   nrows)
        col_tics_start   = gammalib.GFitsTableLongCol('TICS_START',   nrows)
        col_tjd_stop     = gammalib.GFitsTableShortCol('TJD_STOP',    nrows)
        col_tics_stop    = gammalib.GFitsTableLongCol('TICS_STOP',    nrows)
        col_quality      = gammalib.GFitsTableShortCol('QUALITY',     nrows)
        col_oad          = gammalib.GFitsTableStringCol('OAD',        nrows, 60)
        col_superpackets = gammalib.GFitsTableLongCol('SUPERPACKETS', nrows)

        # Set column units
        col_glon_scz.unit('deg')
        col_glat_scz.unit('deg')
        col_tjd_start.unit('days')
        col_tics_start.unit('0.000125 s')
        col_tjd_stop.unit('days')
        col_tics_stop.unit('0.000125 s')

        # Fill columns
        for i, oad_record in enumerate(oad_list):
            col_vp[i]           = oad_record['vp']
            col_object[i]       = oad_record['object']
            col_obs_id[i]       = '%.1f' % oad_record['obs_id']
            col_dsd_id[i]       = oad_record['dsd_id']
            col_glon_scz[i]     = oad_record['glon_scz']
            col_glat_scz[i]     = oad_record['glat_scz']
            col_date_start[i]   = oad_record['date_start']
            col_date_stop[i]    = oad_record['date_stop']
            col_tjd_start[i]    = oad_record['tjd_start']
            col_tics_start[i]   = oad_record['tics_start']
            col_tjd_stop[i]     = oad_record['tjd_stop']
            col_tics_stop[i]    = oad_record['tics_stop']
            col_quality[i]      = oad_record['quality']
            col_oad[i]          = oad_record['filename']
            col_superpackets[i] = oad_record['superpackets']

        # Create FITS table
        table = gammalib.GFitsBinTable()

        # Append columns to table
        table.append(col_vp)
        table.append(col_object)
        table.append(col_obs_id)
        table.append(col_dsd_id)
        table.append(col_glon_scz)
        table.append(col_glat_scz)
        table.append(col_date_start)
        table.append(col_date_stop)
        table.append(col_tjd_start)
        table.append(col_tics_start)
        table.append(col_tjd_stop)
        table.append(col_tics_stop)
        table.append(col_quality)
        table.append(col_oad)
        table.append(col_superpackets)

        # Set extension name
        table.extname('OAD')

        # Return table
        return table

    def _build_hkd_table(self, hkd_list):
        """
        Build HKD FITS table

        Parameters
        ----------
        hkd_list : list of dict
            HKD database

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            HKD FITS table
        """
        # Determine number of rows in table
        nrows = len(hkd_list)

        # Create FITS table columns
        col_vp         = gammalib.GFitsTableStringCol('VP',         nrows, 9)
        col_object     = gammalib.GFitsTableStringCol('OBJECT',     nrows, 20)
        col_obs_id     = gammalib.GFitsTableStringCol('OBS_ID',     nrows, 9)
        col_dsd_id     = gammalib.GFitsTableStringCol('DSD_ID',     nrows, 15)
        col_glon_scz   = gammalib.GFitsTableFloatCol('GLON_SCZ',    nrows)
        col_glat_scz   = gammalib.GFitsTableFloatCol('GLAT_SCZ',    nrows)
        col_date_start = gammalib.GFitsTableStringCol('DATE_START', nrows, 20)
        col_date_stop  = gammalib.GFitsTableStringCol('DATE_STOP',  nrows, 20)
        col_tjd_start  = gammalib.GFitsTableShortCol('TJD_START',   nrows)
        col_tics_start = gammalib.GFitsTableLongCol('TICS_START',   nrows)
        col_tjd_stop   = gammalib.GFitsTableShortCol('TJD_STOP',    nrows)
        col_tics_stop  = gammalib.GFitsTableLongCol('TICS_STOP',    nrows)
        col_quality    = gammalib.GFitsTableShortCol('QUALITY',     nrows)
        col_hkd        = gammalib.GFitsTableStringCol('HKD',        nrows, 60)
        col_records    = gammalib.GFitsTableLongCol('RECORDS',      nrows)

        # Set column units
        col_glon_scz.unit('deg')
        col_glat_scz.unit('deg')
        col_tjd_start.unit('days')
        col_tics_start.unit('0.000125 s')
        col_tjd_stop.unit('days')
        col_tics_stop.unit('0.000125 s')

        # Fill columns
        for i, hkd_record in enumerate(hkd_list):
            col_vp[i]         = hkd_record['vp']
            col_object[i]     = hkd_record['object']
            col_obs_id[i]     = '%.1f' % hkd_record['obs_id']
            col_dsd_id[i]     = hkd_record['dsd_id']
            col_glon_scz[i]   = hkd_record['glon_scz']
            col_glat_scz[i]   = hkd_record['glat_scz']
            col_date_start[i] = hkd_record['date_start']
            col_date_stop[i]  = hkd_record['date_stop']
            col_tjd_start[i]  = hkd_record['tjd_start']
            col_tics_start[i] = hkd_record['tics_start']
            col_tjd_stop[i]   = hkd_record['tjd_stop']
            col_tics_stop[i]  = hkd_record['tics_stop']
            col_quality[i]    = hkd_record['quality']
            col_hkd[i]        = hkd_record['filename']
            col_records[i]    = hkd_record['records']

        # Create FITS table
        table = gammalib.GFitsBinTable()

        # Append columns to table
        table.append(col_vp)
        table.append(col_object)
        table.append(col_obs_id)
        table.append(col_dsd_id)
        table.append(col_glon_scz)
        table.append(col_glat_scz)
        table.append(col_date_start)
        table.append(col_date_stop)
        table.append(col_tjd_start)
        table.append(col_tics_start)
        table.append(col_tjd_stop)
        table.append(col_tics_stop)
        table.append(col_quality)
        table.append(col_hkd)
        table.append(col_records)

        # Set extension name
        table.extname('HKD')

        # Return table
        return table

    def _build_tim_table(self, tim_list):
        """
        Build TIM FITS table

        Parameters
        ----------
        tim_list : list of dict
            TIM database

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            TIM FITS table
        """
        # Determine number of rows in table
        nrows = len(tim_list)

        # Create FITS table columns
        col_vp         = gammalib.GFitsTableStringCol('VP',         nrows, 9)
        col_object     = gammalib.GFitsTableStringCol('OBJECT',     nrows, 20)
        col_obs_id     = gammalib.GFitsTableStringCol('OBS_ID',     nrows, 9)
        col_dsd_id     = gammalib.GFitsTableStringCol('DSD_ID',     nrows, 15)
        col_glon_scz   = gammalib.GFitsTableFloatCol('GLON_SCZ',    nrows)
        col_glat_scz   = gammalib.GFitsTableFloatCol('GLAT_SCZ',    nrows)
        col_date_start = gammalib.GFitsTableStringCol('DATE_START', nrows, 20)
        col_date_stop  = gammalib.GFitsTableStringCol('DATE_STOP',  nrows, 20)
        col_tjd_start  = gammalib.GFitsTableShortCol('TJD_START',   nrows)
        col_tics_start = gammalib.GFitsTableLongCol('TICS_START',   nrows)
        col_tjd_stop   = gammalib.GFitsTableShortCol('TJD_STOP',    nrows)
        col_tics_stop  = gammalib.GFitsTableLongCol('TICS_STOP',    nrows)
        col_telapse    = gammalib.GFitsTableDoubleCol('TELAPSE',    nrows)
        col_ontime     = gammalib.GFitsTableDoubleCol('ONTIME',     nrows)
        col_quality    = gammalib.GFitsTableShortCol('QUALITY',     nrows)
        col_tim        = gammalib.GFitsTableStringCol('TIM',        nrows, 60)
        col_intervals  = gammalib.GFitsTableLongCol('INTERVALS',    nrows)

        # Set column units
        col_glon_scz.unit('deg')
        col_glat_scz.unit('deg')
        col_tjd_start.unit('days')
        col_tics_start.unit('0.000125 s')
        col_tjd_stop.unit('days')
        col_tics_stop.unit('0.000125 s')
        col_telapse.unit('s')
        col_ontime.unit('s')

        # Fill columns
        for i, tim_record in enumerate(tim_list):
            col_vp[i]         = tim_record['vp']
            col_object[i]     = tim_record['object']
            col_obs_id[i]     = '%.1f' % tim_record['obs_id']
            col_dsd_id[i]     = tim_record['dsd_id']
            col_glon_scz[i]   = tim_record['glon_scz']
            col_glat_scz[i]   = tim_record['glat_scz']
            col_date_start[i] = tim_record['date_start']
            col_date_stop[i]  = tim_record['date_stop']
            col_tjd_start[i]  = tim_record['tjd_start']
            col_tics_start[i] = tim_record['tics_start']
            col_tjd_stop[i]   = tim_record['tjd_stop']
            col_tics_stop[i]  = tim_record['tics_stop']
            col_telapse[i]    = tim_record['telapse']
            col_ontime[i]     = tim_record['ontime']
            col_quality[i]    = tim_record['quality']
            col_tim[i]        = tim_record['filename']
            col_intervals[i]  = tim_record['intervals']

        # Create FITS table
        table = gammalib.GFitsBinTable()

        # Append columns to table
        table.append(col_vp)
        table.append(col_object)
        table.append(col_obs_id)
        table.append(col_dsd_id)
        table.append(col_glon_scz)
        table.append(col_glat_scz)
        table.append(col_date_start)
        table.append(col_date_stop)
        table.append(col_tjd_start)
        table.append(col_tics_start)
        table.append(col_tjd_stop)
        table.append(col_tics_stop)
        table.append(col_telapse)
        table.append(col_ontime)
        table.append(col_quality)
        table.append(col_tim)
        table.append(col_intervals)

        # Set extension name
        table.extname('TIM')

        # Return table
        return table

    def _build_xml_table(self, xml_list):
        """
        Build XML FITS table

        Parameters
        ----------
        xml_list : list of dict
            XML database

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            XML FITS table
        """
        # Determine number of rows in table
        nrows = len(xml_list)

        # Create FITS table columns
        col_vp          = gammalib.GFitsTableStringCol('VP',         nrows, 9)
        col_object      = gammalib.GFitsTableStringCol('OBJECT',     nrows, 20)
        col_obs_id      = gammalib.GFitsTableStringCol('OBS_ID',     nrows, 9)
        col_glon_scz    = gammalib.GFitsTableFloatCol('GLON_SCZ',    nrows)
        col_glat_scz    = gammalib.GFitsTableFloatCol('GLAT_SCZ',    nrows)
        col_date_start  = gammalib.GFitsTableStringCol('DATE_START', nrows, 20)
        col_date_stop   = gammalib.GFitsTableStringCol('DATE_STOP',  nrows, 20)
        col_tjd_start   = gammalib.GFitsTableShortCol('TJD_START',   nrows)
        col_tics_start  = gammalib.GFitsTableLongCol('TICS_START',   nrows)
        col_tjd_stop    = gammalib.GFitsTableShortCol('TJD_STOP',    nrows)
        col_tics_stop   = gammalib.GFitsTableLongCol('TICS_STOP',    nrows)
        col_visday      = gammalib.GFitsTableShortCol('VISDAY',      nrows)
        col_vistim      = gammalib.GFitsTableLongCol('VISTIM',       nrows)
        col_vieday      = gammalib.GFitsTableShortCol('VIEDAY',      nrows)
        col_vietim      = gammalib.GFitsTableLongCol('VIETIM',       nrows)
        col_telapse     = gammalib.GFitsTableDoubleCol('TELAPSE',    nrows)
        col_ontime      = gammalib.GFitsTableDoubleCol('ONTIME',     nrows)
        col_evp_verno   = gammalib.GFitsTableShortCol('EVP_VERNO',   nrows)
        col_evp_quality = gammalib.GFitsTableShortCol('EVP_QUALITY', nrows)
        col_tim_quality = gammalib.GFitsTableShortCol('TIM_QUALITY', nrows)
        col_oad_quality = gammalib.GFitsTableShortCol('OAD_QUALITY', nrows, 43)
        col_hkd_quality = gammalib.GFitsTableShortCol('HKD_QUALITY', nrows, 43)
        col_xml         = gammalib.GFitsTableStringCol('XML',        nrows, 60)
        col_evp         = gammalib.GFitsTableStringCol('EVP',        nrows, 60)
        col_tim         = gammalib.GFitsTableStringCol('TIM',        nrows, 60)
        col_oad         = gammalib.GFitsTableStringCol('OAD',        nrows, 60, 43)
        col_hkd         = gammalib.GFitsTableStringCol('HKD',        nrows, 60, 43)
        col_nevents     = gammalib.GFitsTableLongCol('NEVENTS',      nrows)
        col_model       = gammalib.GFitsTableStringCol('MODE',       nrows, 8)

        # Set column units
        col_glon_scz.unit('deg')
        col_glat_scz.unit('deg')
        col_tjd_start.unit('days')
        col_tics_start.unit('0.000125 s')
        col_tjd_stop.unit('days')
        col_tics_stop.unit('0.000125 s')
        col_visday.unit('days')
        col_vistim.unit('0.000125 s')
        col_vieday.unit('days')
        col_vietim.unit('0.000125 s')
        col_telapse.unit('s')
        col_ontime.unit('s')

        # Fill columns
        for i, xml_record in enumerate(xml_list):
            col_vp[i]          = xml_record['vp']
            col_object[i]      = xml_record['object']
            col_obs_id[i]      = '%.1f' % xml_record['obs_id']
            col_glon_scz[i]    = xml_record['glon_scz']
            col_glat_scz[i]    = xml_record['glat_scz']
            col_date_start[i]  = xml_record['date_start']
            col_date_stop[i]   = xml_record['date_stop']
            col_tjd_start[i]   = xml_record['tjd_start']
            col_tics_start[i]  = xml_record['tics_start']
            col_tjd_stop[i]    = xml_record['tjd_stop']
            col_tics_stop[i]   = xml_record['tics_stop']
            col_visday[i]      = xml_record['visday']
            col_vistim[i]      = xml_record['vistim']
            col_vieday[i]      = xml_record['vieday']
            col_vietim[i]      = xml_record['vietim']
            col_telapse[i]     = xml_record['telapse']
            col_ontime[i]      = xml_record['ontime']
            col_evp_verno[i]   = xml_record['evp_verno']
            col_evp_quality[i] = xml_record['evp_quality']
            col_tim_quality[i] = xml_record['tim_quality']
            col_xml[i]         = xml_record['xmlname']
            col_evp[i]         = xml_record['evpname']
            col_tim[i]         = xml_record['timname']
            col_nevents[i]     = xml_record['nevents']
            for k, oad_quality in enumerate(xml_record['oad_quality']):
                col_oad_quality[i,k] = oad_quality
            for k, oadname in enumerate(xml_record['oadnames']):
                col_oad[i,k] = oadname
            for k, hkd_quality in enumerate(xml_record['hkd_quality']):
                col_hkd_quality[i,k] = hkd_quality
            for k, hkdname in enumerate(xml_record['hkdnames']):
                col_hkd[i,k] = hkdname
            mode = 'UNKNOWN'
            if xml_record['mode'] == 'STD':
                mode = 'STANDARD'
            elif xml_record['mode'] == 'S80':
                mode = 'SOLAR80'
            col_model[i] = mode

        # Create FITS table
        table = gammalib.GFitsBinTable()

        # Append columns to table
        table.append(col_vp)
        table.append(col_object)
        table.append(col_obs_id)
        table.append(col_glon_scz)
        table.append(col_glat_scz)
        table.append(col_date_start)
        table.append(col_date_stop)
        table.append(col_tjd_start)
        table.append(col_tics_start)
        table.append(col_tjd_stop)
        table.append(col_tics_stop)
        table.append(col_visday)
        table.append(col_vistim)
        table.append(col_vieday)
        table.append(col_vietim)
        table.append(col_telapse)
        table.append(col_ontime)
        table.append(col_nevents)
        table.append(col_evp_verno)
        table.append(col_evp_quality)
        table.append(col_tim_quality)
        table.append(col_oad_quality)
        table.append(col_hkd_quality)
        table.append(col_xml)
        table.append(col_evp)
        table.append(col_tim)
        table.append(col_oad)
        table.append(col_hkd)

        # Set extension name
        table.extname('XML')

        # Return table
        return table


    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Optionally download COMPTEL HEASARC archive
        if self['download'].boolean():
            self._download()

        # Log header
        self._log_header1(gammalib.TERSE, 'Get list of EVP files')

        # Get EVP list
        self._evp_list = self._get_evp_list()

        # Log header
        self._log_header1(gammalib.TERSE, 'Get list of OAD files')

        # Get OAD list
        self._oad_list = self._get_oad_list()

        # Log header
        self._log_header1(gammalib.TERSE, 'Get list of HKD files')

        # Get HKD list
        self._hkd_list = self._get_hkd_list()

        # Log header
        self._log_header1(gammalib.TERSE, 'Get list of TIM files')

        # Get TIM list
        self._tim_list = self._get_tim_list()

        # Log header
        self._log_header1(gammalib.TERSE, 'Generate list of XML files')

        # Get XML list
        xml_list = self._get_xml_list(self._evp_list, self._tim_list, self._oad_list, self._hkd_list)

        # Log header
        self._log_header1(gammalib.TERSE, 'Compare with COMPASS database')

        # Check XML list
        self._xml_list = self._check_xml_list(xml_list)

        # Return
        return

    def save(self):
        """ 
        Save database file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save database')

        # Create database directory structure
        self._create_dbase_dirs()

        # Generate TIM and XML files
        for xml in self._xml_list:
            self._build_tim_file(xml)
            self._build_xml_file(xml)

        # Check database
        self._check_database(self._xml_list)

        # Build database file
        self._build_database_file(self._evp_list, self._oad_list,
                                  self._hkd_list, self._tim_list,
                                  self._xml_list)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comgendb(sys.argv)

    # Execute application
    app.execute()
