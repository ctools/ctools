#!/usr/bin/env python
# ==========================================================================
# Generation of an caldb entry from general IACT observation.
#
# Copyright (C) 2015-2016 Michael Mayer
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
import os
import datetime


# ================= #
# csobs2caldb class #
# ================= #
class csobs2caldb(ctools.cscript):
    """
    Creates a calibration database entry for an observation.
    
    The creation of a calibration database entry is useful for performing
    simulations for current Imaging Air Cherenkov Telescopes (IACTs).
    The class takes an observation definition XML file as input and uses
    one observation to create a calibration database entry. By default
    the first observation will be used, but the user can specify the
    index of any observation using the hidden "index" parameter.
    """
    
    # Constructors and destructors
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name and version
        self._name    = "csobs2caldb"
        self._version = "1.1.0"

        # Initialise members
        self._obs         = []     # Observation container
        self._observation = None
        self._mission     = "CTA"
        self._instrument  = ""
        self._rspname     = ""
        self._outfile     = ""
        self._base_dir    = ""
        self._cal_dir     = ""
        
        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self._name, self._version)
        elif len(argv) == 1:
            ctools.cscript.__init__(self, self._name, self._version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self._log_header()
        self._log.date(True)

        # Return
        return

    def __del__(self):
        """
        Destructor.
        """
        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Load observation container
        self._obs = gammalib.GObservations(self["inobs"].filename())
        
        # Get index of observation to be used (0=default)
        inx = self["index"].integer()
        
        # Get observation
        self._observation = self._obs[inx]
        
        # Get instrument name
        self._instrument = self._observation.instrument()
        
        # Get config and file name
        self._rspname = self["rspname"].string()
        self._outfile = self["outirfs"].filename()
        
        # Make sure we have a cta observation
        if not self._observation.classname() == "GCTAObservation":
            raise RuntimeError("Input observation not of type "+
                               "\"GCTAObservation\".")
        
        # Make sure we have an irf response
        if not self._observation.response().classname() == "GCTAResponseIrf":
            raise RuntimeError("Response of input observation not of "+
                               "type \"GCTAResponseIrf\".")

        # Get other parameters
        self._clobber = self["clobber"].boolean()

        # Return
        return

    def _make_irf_file(self):
        """
        Creates an IRF FITS file.
        """
        # Write header into logger
        if self._logTerse():
            self._log("\n")
            self._log.header2("Creating IRF file")

        # Get response for the observation
        rsp = self._observation.response()
        
        # Extract response file names
        fname_aeff  = rsp.aeff().filename()
        fname_psf   = rsp.psf().filename()
        fname_edisp = rsp.edisp().filename()
        fname_bkg   = rsp.background().filename()

        # Log filename instances
        if self._logTerse():
            self._log(gammalib.parformat("Effective area"))
            self._log(repr(fname_aeff))
            self._log("\n")
            self._log(gammalib.parformat("Point spread function"))
            self._log(fname_psf.url())
            self._log("\n")
            self._log(gammalib.parformat("Energy dispersion"))
            self._log(fname_edisp.url())
            self._log("\n")
            self._log(gammalib.parformat("Background rate"))
            self._log(fname_bkg.url())
            self._log("\n")
    
        # Open FITS files of response components
        fits_aeff  = gammalib.GFits(fname_aeff)
        fits_psf   = gammalib.GFits(fname_psf)
        fits_edisp = gammalib.GFits(fname_edisp)
        fits_bkg   = gammalib.GFits(fname_bkg)
        
        # Create empty FITS file
        fits = gammalib.GFits()
        
        # Append IRF component to FITS file
        fits.append(fits_aeff[fname_aeff.extname("EFFECTIVE AREA")])
        fits.append(fits_psf[fname_psf.extname("POINT SPREAD FUNCTION")])
        fits.append(fits_edisp[fname_edisp.extname("ENERGY DISPERSION")])
        fits.append(fits_bkg[fname_bkg.extname("BACKGROUND")])
        
        # Return fits file
        return fits

    def _make_caldb_index(self):
        """
        Creates an IRF index FITS file.
        """
        # Write header into logger
        if self._logTerse():
            self._log("\n")
            self._log.header2("Creating database index file")

        # Open calibration database index
        indx_file = self._base_dir+"/caldb.indx"
        
        # Open index file (or create one if it does not exist)
        cif = gammalib.GFits(indx_file, True)
        
        # Retrieve "CIF" table
        if cif.contains("CIF"):
            table = cif.table("CIF")

        # ... or create binary table if no "CIF" table exists yet
        else:       

            # Create empty binary table
            bintable = gammalib.GFitsBinTable()
            bintable.extname("CIF")
            bintable.card("CIFVERSN", "1992a", "Version of CIF format")     

            # Append table to FITS file and recover a reference to the
            # appended table
            table = cif.append(bintable)
            
            # Attach columns. Reference: CAL/GEN/92-008
            table.append(gammalib.GFitsTableStringCol("TELESCOP", 0, 10))
            table.append(gammalib.GFitsTableStringCol("INSTRUME", 0, 10))
            table.append(gammalib.GFitsTableStringCol("DETNAM", 0, 20))
            table.append(gammalib.GFitsTableStringCol("FILTER", 0, 10))
            table.append(gammalib.GFitsTableStringCol("CAL_DEV", 0, 20))
            table.append(gammalib.GFitsTableStringCol("CAL_DIR", 0, 70))
            table.append(gammalib.GFitsTableStringCol("CAL_FILE", 0, 40))
            table.append(gammalib.GFitsTableStringCol("CAL_CLAS", 0, 3))
            table.append(gammalib.GFitsTableStringCol("CAL_DTYP", 0, 4))
            table.append(gammalib.GFitsTableStringCol("CAL_CNAM", 0, 20))
            table.append(gammalib.GFitsTableStringCol("CAL_CBD", 0, 70, 9))
            table.append(gammalib.GFitsTableShortCol("CAL_XNO", 0))
            table.append(gammalib.GFitsTableStringCol("CAL_VSD", 0, 10))
            table.append(gammalib.GFitsTableStringCol("CAL_VST", 0, 8))
            table.append(gammalib.GFitsTableDoubleCol("REF_TIME", 0))
            table.append(gammalib.GFitsTableShortCol("CAL_QUAL", 0))
            table.append(gammalib.GFitsTableStringCol("CAL_DATE", 0, 8))
            table.append(gammalib.GFitsTableStringCol("CAL_DESC", 0, 70))

        # Check if output config already exist
        has_config = False
        row_index = -1
        for caldir in table["CAL_DIR"]:
            row_index += 1
            if caldir == self._cal_dir:
                has_config = True
                break
         
        # Create columns if not available   
        if not has_config:           
            # Append 4 rows to CIF extension
            table.append_rows(4)
            row_index = table.nrows()
        else:
            row_index += 4

        # Add generic information for these 4 rows
        for i in range(4):

            # Set row number
            row = i + row_index-4

            # Set date
            now = str(datetime.datetime.now())
            today, time = now.split()
            
            # Set element
            table["TELESCOP"][row] = self._mission
            table["INSTRUME"][row] = self._instrument
            table["DETNAM"][row]   = "NONE"
            table["FILTER"][row]   = "NONE"
            table["CAL_DEV"][row]  = "ONLINE"
            table["CAL_CLAS"][row] = "BCF"
            table["CAL_DTYP"][row] = "DATA"
            table["CAL_VSD"][row]  = today
            table["CAL_VST"][row]  = time.split(".")[0]
            table["REF_TIME"][row] = 51544.0
            table["CAL_QUAL"][row] = 0
            table["CAL_CBD"][row]  = "NAME("+self._rspname+")"
            table["CAL_DATE"][row] = today.replace("-","/")[2:]
            table["CAL_DIR"][row]  = self._cal_dir
            table["CAL_FILE"][row] = self._outfile

        # Add effective area information
        row = row_index-4
        table["CAL_CNAM"][row] = "EFF_AREA"
        table["CAL_DESC"][row] = self._instrument+" effective area"

        # Add point spread function information
        row = row_index-3
        table["CAL_CNAM"][row] = "RPSF"
        table["CAL_DESC"][row] = self._instrument+" point spread function"
        
        # Add energy dispersion information
        row = row_index-2
        table["CAL_CNAM"][row] = "EDISP"
        table["CAL_DESC"][row] = self._instrument+" energy dispersion"

        # Add background information
        row = row_index-1
        table["CAL_CNAM"][row] = "BGD"
        table["CAL_DESC"][row] = self._instrument+" background"
        
        # Return CIF FITS file
        return cif
        
    def _make_dirs(self):
        """
        Make CALDB directories.
        """
        # Write header into logger
        if self._logTerse():
            self._log("\n")
            self._log.header2("Creating directory structure")

        # Create calibration database        
        caldb = gammalib.GCaldb()
        
        # Set calibration directory 
        self._cal_dir  = "data"
        self._cal_dir += "/"+self._mission.lower()
        self._cal_dir += "/"+self._instrument.lower()
        self._cal_dir += "/bcf/"+self._rspname
        
        # Set absolute path
        self._base_dir = caldb.rootdir() +"/data"
        self._base_dir += "/"+self._mission.lower()
        self._base_dir += "/"+self._instrument.lower()
        
        # Set directory for irf file
        self._rsp_dir = caldb.rootdir() + "/" + self._cal_dir
        
        # Create directories and log information
        if not os.path.isdir(self._rsp_dir):
            if self._logExplicit():
                self._log(gammalib.parformat("Directory"))
                self._log(self._rsp_dir)
                self._log("\n")
            os.makedirs(self._rsp_dir)
        else:
            if self._logExplicit():
                self._log(gammalib.parformat("Directory (existing)"))
                self._log(self._rsp_dir)
                self._log("\n")

        # Return
        return


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
        
        # Write input parameters and header into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n\n")
            self._log.header1("Creating CALDB entry")
        
        # Create directory structure
        self._make_dirs()
         
        # Create/update calibration database   
        self._caldb_inx = self._make_caldb_index()
        
        # Create response file
        self._irf_fits = self._make_irf_file()
                
        # Return
        return

    def save(self):
        """
        Save calibration database FITS file.
        """
        # Write header into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save calibration database")

        # Set response filename
        filename = self._rsp_dir + "/" + self._outfile
        
        # Write filenames into logger
        if self._logTerse():
            self._log(gammalib.parformat("CALDB index file"))
            self._log(self._caldb_inx.filename().url())
            self._log("\n")
            self._log(gammalib.parformat("Response file"))
            self._log(filename)
            self._log("\n")

        # Save caldb index file
        self._caldb_inx.save(self._clobber)

        # Save response file
        self._irf_fits.saveto(filename, self._clobber)
        
        # Return
        return
        
    def execute(self):
        """
        Execute the script.
        """
        # Open logfile
        self._logFileOpen()

        # Run the script
        self.run()

        # Save to disk
        self.save()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csobs2caldb(sys.argv)
    
    # Execute application
    app.execute()
