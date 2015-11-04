#! /usr/bin/env python
# ==========================================================================
# Generation of an observation definition file.
#
# Copyright (C) 2015 Juergen Knoedlseder
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
    The csobsdef class generates an observation definition file from a
    pointing list. The pointing list is a comma-separated value ASCII
    file with header keywords in the first row followed by a list of
    pointings (one pointing per row). The following header keywords
    are supported (case sensitive, column order irrelevant):

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

    Here some usage examples:

    csobsdef
      Creates minimal observation definition file.

    csobsdef emin=0.1 emax=100.0
      Creates observation definition file with an energy range 100 GeV - 100 TeV.

    csobsdef rad=5
      Creates observation definition file with a ROI radius of 5 deg.

    csobsdef caldb=dummy irf=cta_dummy_irf
      Creates observation definition file using the "cta_dummy_irf" IRF in the
      "dummy" calibration database.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name and version
        self.name    = "csobsdef"
        self.version = "1.0.0"

        # Initialise class members
        self.obs        = gammalib.GObservations()
        self.inpnt      = ""
        self.tmin       = 0.0
        self.outfile    = ""
        self.chatter    = 2
        self.clobber    = True
        self.debug      = False
        self.read_ahead = False

        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self.name, self.version)
        elif len(argv) ==1:
            ctools.cscript.__init__(self, self.name, self.version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self.log_header()
        self.log.date(True)

        # Return
        return

    def __del__(self):
        """
        Destructor.
        """
        # Return
        return

    def parfile(self):
        """
        Check if parfile exists. If parfile does not exist then create a
        default parfile. This kluge avoids shipping the cscript with a parfile.
        """
        # Set parfile name
        parfile = self.name+".par"

        # Test parfile existence, and create parfile if it does not exist
        try:
            pars = gammalib.GApplicationPars(parfile)
        except:
            # Signal that parfile was not found
            print("Parfile \""+parfile+"\" not found. Create default parfile.")

            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("inpnt","f","a","NONE","","","Pointing definition file"))
            pars.append(gammalib.GApplicationPar("outobs","f","a","obs.xml","","","Output observation definition file"))
            pars.append(gammalib.GApplicationPar("duration","r","a","1800.0","","","Pointing duration (seconds)"))
            pars.append(gammalib.GApplicationPar("caldb","s","h","","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","h","","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("emin","r","h","UNDEF","","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","h","UNDEF","","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("rad","r","h","UNDEF","","","Radius of ROI (degrees)"))
            pars.append(gammalib.GApplicationPar("deadc","r","h","0.95","","","Deadtime correction factor"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csobsdef.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile.
        """
        # Read parameters    
        self.inpnt = self["inpnt"].filename()

        # Read ahead parameters
        if self.read_ahead:
            self.outfile = self["outobs"].filename()

        # Set some fixed parameters
        self.chatter = self["chatter"].integer()
        self.clobber = self["clobber"].boolean()
        self.debug   = self["debug"].boolean()

        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Set read ahead flag
        self.read_ahead = True

        # Run the script
        self.run()

        # Save observation definition file
        self.obs.save(self.outfile)

        # Return
        return

    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self.logDebug():
            self.log.cout(True)

        # Get parameters
        self.get_parameters()

        # Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")

        # Load pointing definition file
        pntdef = gammalib.GCsv(self.inpnt, ",")
        ncols  = pntdef.ncols()
        npnt   = pntdef.nrows()-1

        # Throw an exception is there is no header information
        if pntdef.nrows() < 1:
            raise RuntimeError("No header found in pointing definition file.")

        # Clear observation container
        self.obs.clear()
        self.id = 1

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
                id       = "%6.6d" % self.id
                self.id += 1
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
                raise gammalib.GException.invalid_value(self.name,
                      "No (ra,dec) or (lon,lat) columns found in pointing"
                      " definition file.")
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
                obs = self.__set_response(obs, caldb, irf)

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
            tmin      = self.tmin
            tmax      = self.tmin + duration
            gti       = gammalib.GGti(self.time_reference())
            tstart    = gammalib.GTime(tmin, self.time_reference())
            tstop     = gammalib.GTime(tmax, self.time_reference())
            self.tmin = tmax
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

            # Append observation
            self.obs.append(obs)

        # Return
        return

    def __set_response(self, obs, caldb, irf):
        """
        Set response for an observation. We create an XML element for that
        so that the response XML writer will write the database and response
        name into the observation definition file.
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


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates observation definition file.
    """
    # Create instance of application
    app = csobsdef(sys.argv)

    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
