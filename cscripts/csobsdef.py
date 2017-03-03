#! /usr/bin/env python
# ==========================================================================
# Generation of an observation definition file
#
# Copyright (C) 2015-2017 Juergen Knoedlseder
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


# ============== #
# csobsdef class #
# ============== #
class csobsdef(ctools.cscript):
    """
    Creates an observation definition file from a pointing list
    
    The pointing list is a comma-separated value ASCII file with header
    keywords in the first row followed by a list of pointings (one
    pointing per row). The following header keywords are supported (case
    sensitive, column order irrelevant):

    name     - Observation name string
    id       - Unique observation identifier string
    ra       - Right Ascension of pointing (deg)
    dec      - Declination of pointing (deg)
    lon      - Galactic longitude of pointing (deg)
    lat      - Galactic latitude of pointing (deg)
    tmin     - Start of pointing (seconds)
    duration - Duration of pointing (seconds)
    emin     - Lower energy limit (TeV)
    emax     - Upper energy limit (TeV)
    rad      - Radius of region of interest (deg)
    deadc    - Deadtime correction factor [0-1]
    caldb    - Calibration database
    irf      - Response function name

    Only the pairs (ra,dec) or (lon,lat) are mandatory header keywords.
    All other keywords are optional and can be specified when calling
    csobsdef as user parameters. The only exception is the "duration"
    keyword that will automatically be queried.

    Examples:

        ./csobsdef
            Creates minimal observation definition file.

        ./csobsdef emin=0.1 emax=100.0
            Creates observation definition file with an energy range
            100 GeV - 100 TeV.

        ./csobsdef rad=5
            Creates observation definition file with a ROI radius of
            5 degrees.

        ./csobsdef caldb=prod2 irf=South_50h
            Creates observation definition file using the "South_50h"
            IRF in the "prod2" calibration database.
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
        """
        # Set name and version
        self._name    = 'csobsdef'
        self._version = '1.3.0'

        # Initialise class members
        self._obs    = gammalib.GObservations()
        self._pntdef = gammalib.GCsv()
        self._tmin   = 0.0

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query input filename if necessary
        if self._pntdef.size() == 0:
            self['inpnt'].filename()

        # Read ahead parameters
        if self._read_ahead():
            self['outobs'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _set_response(self, obs, caldb, irf):
        """
        Set response for an observation
        
        The method creates an XML element so that that the response XML
        writer will write the database and response name into the
        observation definition file.

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation`
            CTA observation
        caldb : str
            Calibration database
        irf : str
            Instrument response function

        Returns
        obs : `~gammalib.GCTAObservation`
            CTA observation with response attached
        -------
        """
        # Create XML element
        xml = gammalib.GXmlElement()

        # Append parameter
        parameter = 'parameter name="Calibration" database="'+caldb+\
                    '" response="'+irf+'"'
        xml.append(gammalib.GXmlElement(parameter))

        # Create CTA response
        response = gammalib.GCTAResponseIrf()
        response.read(xml)

        # Attach response to observation
        obs.response(response)

        # Return observation
        return obs

 
    # Public methods
    def run(self):
        """
        Run the script

        Raises
        ------
        RuntimeError
            Invalid pointing definition file format
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Write header into logger
        self._log_header1(gammalib.TERSE,
                          'Creating observation definition XML file')

        # Load pointing definition file if it is not already set. Extract
        # the number of columns and pointings
        if self._pntdef.size() == 0:
            self._pntdef = gammalib.GCsv(self['inpnt'].filename(), ',')
        ncols = self._pntdef.ncols()
        npnt  = self._pntdef.nrows()-1

        # Raise an exception if there is no header information
        if self._pntdef.nrows() < 1:
            raise RuntimeError('No header found in pointing definition file.')

        # Clear observation container
        self._obs.clear()

        # Initialise observation identifier counter
        identifier = 1

        # Extract header columns from pointing definition file and put them
        # into a list
        header = []
        for col in range(ncols):
            header.append(self._pntdef[0,col])

        # Loop over all pointings
        for pnt in range(npnt):

            # Set pointing definition CSV file row index
            row = pnt + 1

            # Create empty CTA observation
            obs = gammalib.GCTAObservation()

            # Set observation name. If no observation name was given then
            # use "None".
            if 'name' in header:
                name = self._pntdef[row, header.index('name')]
            else:
                name = self['name'].string()
            obs.name(name)

            # Set observation identifier. If no observation identified was
            # given the use the internal counter.
            if 'id' in header:
                obsid = self._pntdef[row, header.index('id')]
            else:
                obsid = '%6.6d' % identifier
                identifier += 1
            obs.id(obsid)

            # Set pointing. Either use "ra" and "dec" or "lon" and "lat".
            # If none of these pairs are given then raise an exception.
            if 'ra' in header and 'dec' in header:
                ra     = float(self._pntdef[row, header.index('ra')])
                dec    = float(self._pntdef[row, header.index('dec')])
                pntdir = gammalib.GSkyDir()
                pntdir.radec_deg(ra,dec)
            elif 'lon' in header and 'lat' in header:
                lon    = float(self._pntdef[row, header.index('lon')])
                lat    = float(self._pntdef[row, header.index('lat')])
                pntdir = gammalib.GSkyDir()
                pntdir.lb_deg(lon,lat)
            else:
                raise RuntimeError('No (ra,dec) or (lon,lat) columns '
                                   'found in pointing definition file.')
            obs.pointing(gammalib.GCTAPointing(pntdir))

            # Set response function. If no "caldb" or "irf" information is
            # provided then use the user parameter values.
            if 'caldb' in header:
                caldb = self._pntdef[row, header.index('caldb')]
            else:
                caldb = self['caldb'].string()
            if 'irf' in header:
                irf = self._pntdef[row, header.index('irf')]
            else:
                irf = self['irf'].string()
            if caldb != '' and irf != '':
                obs = self._set_response(obs, caldb, irf)

            # Set deadtime correction factor. If no information is provided
            # then use the user parameter value "deadc".
            if 'deadc' in header:
                deadc = float(self._pntdef[row, header.index('deadc')])
            else:
                deadc = self['deadc'].real()
            obs.deadc(deadc)

            # Set Good Time Interval. If no information is provided then use
            # the user parameter values "tmin" and "duration".
            if 'tmin' in header:
                self._tmin = float(self._pntdef[row, header.index('tmin')])
            if 'duration' in header:
                duration = float(self._pntdef[row, header.index('duration')])
            else:
                duration = self['duration'].real()
            tmin       = self._tmin
            tmax       = self._tmin + duration
            gti        = gammalib.GGti(self._time_reference())
            tstart     = gammalib.GTime(tmin, self._time_reference())
            tstop      = gammalib.GTime(tmax, self._time_reference())
            self._tmin = tmax
            gti.append(tstart, tstop)
            obs.ontime(gti.ontime())
            obs.livetime(gti.ontime()*deadc)

            # Set Energy Boundaries. If no "emin" or "emax" information is
            # provided then use the user parameter values in case they are
            # valid.
            has_emin = False
            has_emax = False
            if 'emin' in header:
                emin     = float(self._pntdef[row, header.index('emin')])
                has_emin = True
            else:
                if self['emin'].is_valid():
                    emin     = self['emin'].real()
                    has_emin = True
            if 'emax' in header:
                emax     = float(self._pntdef[row, header.index('emax')])
                has_emax = True
            else:
                if self['emax'].is_valid():
                    emax     = self['emax'].real()
                    has_emax = True
            has_ebounds = has_emin and has_emax
            if has_ebounds:
                ebounds = gammalib.GEbounds(gammalib.GEnergy(emin, 'TeV'),
                                            gammalib.GEnergy(emax, 'TeV'))

            # Set ROI. If no ROI radius is provided then use the user
            # parameters "rad".
            has_roi = False
            if 'rad' in header:
                rad     = float(self._pntdef[row, header.index('rad')])
                has_roi = True
            else:
                if self['rad'].is_valid():
                    rad     = self['rad'].real()
                    has_roi = True
            if has_roi:
                roi = gammalib.GCTARoi(gammalib.GCTAInstDir(pntdir), rad)

            # Create an empty event list
            event_list = gammalib.GCTAEventList()
            event_list.gti(gti)

            # If available, set the energy boundaries and the ROI
            if has_ebounds:
                event_list.ebounds(ebounds)
            if has_roi:
                event_list.roi(roi)

            # Attach event list to CTA observation
            obs.events(event_list)

            # Write observation into logger
            name  = obs.instrument()+' observation'
            value = 'Name="%s" ID="%s"' % (obs.name(), obs.id())
            self._log_value(gammalib.NORMAL, name, value)
            self._log_string(gammalib.EXPLICIT, str(obs)+'\n')

            # Append observation
            self._obs.append(obs)

        # Return
        return

    def save(self):
        """
        Save observation definition XML file.
        """
        # Write header and filename into logger
        self._log_header1(gammalib.TERSE, 'Save observation definition XML file')

        # Get output filename in case it was not read ahead
        outobs = self['outobs'].filename()

        # Check if observation definition XML file is valid
        if outobs.url() != 'NONE':

            # Log filename
            self._log_value(gammalib.NORMAL, 'Observation XML file', outobs.url())

            # Save observation definition XML file
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

    def obs(self):
        """
        Returns observation container
        -----
        Returns
        obs : `~gammalib.GObservations`
            Observation container
        """
        # Return container
        return self._obs

    def pntdef(self, csv):
        """
        Set pointing definition from a CSV table

        Parameters
        ----------
        csv : `~gammalib.GCsv`
            Comma-separated values table
        """
        # Set pointing definition
        self._pntdef = csv

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csobsdef(sys.argv)

    # Execute application
    app.execute()
