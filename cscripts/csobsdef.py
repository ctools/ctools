#! /usr/bin/env python
# ==========================================================================
# Generation of an observation definition file.
#
# Copyright (C) 2015-2016 Juergen Knoedlseder
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
import sys


# ============== #
# csobsdef class #
# ============== #
class csobsdef(ctools.cscript):
    """
    Creates an observation definition file from a pointing list.
    
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

        ./csobsdef caldb=dummy irf=cta_dummy_irf
            Creates observation definition file using the "cta_dummy_irf"
            IRF in the "dummy" calibration database.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name and version
        self._name    = "csobsdef"
        self._version = "1.1.0"

        # Initialise class members
        self._obs        = gammalib.GObservations()
        self._inpnt      = ""
        self._tmin       = 0.0
        self._chatter    = 2
        self._clobber    = True
        self._debug      = False

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile.
        """
        # Read parameters    
        self._inpnt = self["inpnt"].filename()

        # Read ahead parameters
        if self._read_ahead():
            self["outobs"].filename()

        # Set some fixed parameters
        self._chatter = self["chatter"].integer()
        self._clobber = self["clobber"].boolean()
        self._debug   = self["debug"].boolean()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Return
        return

    def _set_response(self, obs, caldb, irf):
        """
        Set response for an observation.
        
        The method creates an XML element so that that the response XML
        writer will write the database and response name into the
        observation definition file.
        """
        # Create XML element
        xml = gammalib.GXmlElement()

        # Append parameter
        parameter = "parameter name=\"Calibration\" database=\""+\
                    caldb+"\" response=\""+irf+"\""
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
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Write header into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Creating observation definition XML file")

        # Load pointing definition file
        pntdef = gammalib.GCsv(self._inpnt, ",")
        ncols  = pntdef.ncols()
        npnt   = pntdef.nrows()-1

        # Throw an exception is there is no header information
        if pntdef.nrows() < 1:
            raise RuntimeError("No header found in pointing definition file.")

        # Clear observation container
        self._obs.clear()
        self._id = 1

        # Extract header from pointing definition file
        header = []
        for col in range(ncols):
            header.append(pntdef[0,col])

        # Loop over all pointings
        for pnt in range(npnt):

            # Set row index
            row = pnt + 1

            # Create CTA observation
            obs = gammalib.GCTAObservation()

            # Set observation name
            if "name" in header:
                name = pntdef[row, header.index("name")]
            else:
                name = "None"
            obs.name(name)

            # Set identifier
            if "id" in header:
                id = pntdef[row, header.index("id")]
            else:
                id        = "%6.6d" % self._id
                self._id += 1
            obs.id(id)

            # Set pointing
            if "ra" in header and "dec" in header:
                ra     = float(pntdef[row, header.index("ra")])
                dec    = float(pntdef[row, header.index("dec")])
                pntdir = gammalib.GSkyDir()
                pntdir.radec_deg(ra,dec)
            elif "lon" in header and "lat" in header:
                lon    = float(pntdef[row, header.index("lon")])
                lat    = float(pntdef[row, header.index("lat")])
                pntdir = gammalib.GSkyDir()
                pntdir.lb_deg(lon,lat)
            else:
                raise RuntimeError("No (ra,dec) or (lon,lat) columns "
                                   "found in pointing definition file.")
            obs.pointing(gammalib.GCTAPointing(pntdir))

            # Set response function
            if "caldb" in header:
                caldb = pntdef[row, header.index("caldb")]
            else:
                caldb = self["caldb"].string()
            if "irf" in header:
                irf = pntdef[row, header.index("irf")]
            else:
                irf = self["irf"].string()
            if caldb != "" and irf != "":
                obs = self._set_response(obs, caldb, irf)

            # Set deadtime correction factor
            if "deadc" in header:
                deadc = float(pntdef[row, header.index("deadc")])
            else:
                deadc = self["deadc"].real()
            obs.deadc(deadc)

            # Set Good Time Interval
            if "duration" in header:
                duration = float(pntdef[row, header.index("duration")])
            else:
                duration = self["duration"].real()
            tmin       = self._tmin
            tmax       = self._tmin + duration
            gti        = gammalib.GGti(self._time_reference())
            tstart     = gammalib.GTime(tmin, self._time_reference())
            tstop      = gammalib.GTime(tmax, self._time_reference())
            self._tmin = tmax
            gti.append(tstart, tstop)
            obs.ontime(gti.ontime())
            obs.livetime(gti.ontime()*deadc)

            # Set Energy Boundaries
            has_emin = False
            has_emax = False
            if "emin" in header:
                emin     = float(pntdef[row, header.index("emin")])
                has_emin = True
            else:
                if self["emin"].is_valid():
                    emin     = self["emin"].real()
                    has_emin = True
            if "emax" in header:
                emax     = float(pntdef[row, header.index("emax")])
                has_emax = True
            else:
                if self["emax"].is_valid():
                    emax     = self["emax"].real()
                    has_emax = True
            has_ebounds = has_emin and has_emax
            if has_ebounds:
                ebounds = gammalib.GEbounds(gammalib.GEnergy(emin, "TeV"),
                                            gammalib.GEnergy(emax, "TeV"))

            # Set ROI
            has_roi = False
            if "rad" in header:
                rad     = float(pntdef[row, header.index("rad")])
                has_roi = True
            else:
                if self["rad"].is_valid():
                    rad     = self["rad"].real()
                    has_roi = True
            if has_roi:
                roi = gammalib.GCTARoi(gammalib.GCTAInstDir(pntdir), rad)

            # Create an empty event list
            list = gammalib.GCTAEventList()
            list.gti(gti)

            # Set optional information
            if has_ebounds:
                list.ebounds(ebounds)
            if has_roi:
                list.roi(roi)

            # Attach event list to CTA observation
            obs.events(list)

            # Write observation into logger
            if self._logExplicit():
                self._log(str(obs))
                self._log("\n")
            elif self._logTerse():
                self._log(gammalib.parformat(obs.instrument()+" observation"))
                self._log("Name=\""+obs.name()+"\" ")
                self._log("ID=\""+obs.id()+"\"\n")

            # Append observation
            self._obs.append(obs)

        # Return
        return

    def save(self):
        """
        Save observation definition XML file.
        """
        # Write header and filename into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save observation definition XML file")

        # Get output filename in case it was not read ahead
        outobs = self["outobs"].filename()

        # Check if observation definition XML file is valid
        if outobs.url() != "NONE":      

            # Log filename
            if self._logTerse():
                self._log(gammalib.parformat("Observation XML file"))
                self._log(outobs.url())
                self._log("\n")

            # Save observation definition XML file
            self._obs.save(outobs)
        
        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Open logfile
        self._logFileOpen()

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
    app = csobsdef(sys.argv)

    # Execute application
    app.execute()
