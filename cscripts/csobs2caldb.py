#!/usr/bin/env python
# ==========================================================================
# Generation of an caldb entry from general IACT observation.
#
# Copyright (C) 2015 Michael Mayer
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
    This class implements the creation of a caldb entry for a particular
    observation, which might be helpful for running simulations for current
    IACTs. It takes an observation definition file as input and uses the
    first observation to create a caldb entry which can be used for
    simulations and sensitivity studies later on.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "csobs2caldb"
        self.version = "1.1.0"

        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self.name, self.version)
        elif len(argv) == 1:
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
        default parfile. This kluge avoids shipping the cscript with a
        parfile.
        """

        # Set parfile name
        parfile = self.name+".par"

        try:
            pars = gammalib.GApplicationPars(parfile)
        except:
            # Signal that parfile was not found
            print("Parfile \""+parfile+"\" not found. Create default parfile.")

            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("inobs","f","a","obs.xml","","","Input observation definition file"))
            pars.append(gammalib.GApplicationPar("rspname","s","a","NONE","","","Response output name (e.g. Zenith50)"))            
            pars.append(gammalib.GApplicationPar("outirfs","f","h","irf_file.fits","","","Output IRF file name"))            
            pars.append(gammalib.GApplicationPar("index","i","h","0","","","Index of observation to be used from XML file"))            
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csobs2caldb.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters     
        self.inobs = self["inobs"].filename()
        self.obs   = gammalib.GObservations(self.inobs)
        
        # Get index of observation to be used (0=default)
        inx = self["index"].integer()
        
        # Get observation
        self.observation = self.obs[inx]
        
        # Set instrument mission to cta
        # This is to be compliant with standard gammalib caldb
        self.m_mission = "CTA"
        
        # Get instrument name
        self.m_instrument = self.observation.instrument()
        
        # Get config and file name
        self.m_rspname = self["rspname"].string()
        self.m_outfile = self["outirfs"].filename()
        
        # Make sure we have a cta observation
        if not self.observation.classname() == "GCTAObservation":
            raise gammalib.Exception.runtime_error("Input observation not of type \"GCTAObservation\"")
        
        # Make sure we have an irf response
        if not self.observation.response().classname() == "GCTAResponseIrf":
            raise gammalib.Exception.runtime_error("Response of input observation not of type \"GCTAResponseIrf\"")

        # Get other parameters
        self.m_clobber = self["clobber"].boolean()

        # Return
        return

    def make_irf_file(self):
        
        # retrieve response component
        rsp = self.observation.response()
    
        # Initialise FITS file
        fits = gammalib.GFits()
    
        # Open FITS files of response components
        fits_aeff  = gammalib.GFits(rsp.aeff().filename())
        fits_psf   = gammalib.GFits(rsp.psf().filename())
        fits_edisp = gammalib.GFits(rsp.edisp().filename())
        fits_bkg   = gammalib.GFits(rsp.background().filename())
        
        # Bundle IRFs into one file
        fits.append(fits_aeff["EFFECTIVE AREA"])
        fits.append(fits_psf["POINT SPREAD FUNCTION"])
        fits.append(fits_edisp["ENERGY DISPERSION"])
        fits.append(fits_bkg["BACKGROUND"])
        
        # Return fits file
        return fits

    def make_caldb_index(self):
        
        # Open calibration database index
        indx_file = self.base_dir+"/caldb.indx"
        
        # Check if indx file exists (create if not existent)
        if os.path.isfile(indx_file):
            cif = gammalib.GFits(indx_file)
        else:
            cif = gammalib.GFits(indx_file, True)
        
        # Retrieve "CIF" table
        if cif.contains("CIF"):
            table = cif.table("CIF")
        else:       
            # Create binary table if not available
            bintable = gammalib.GFitsBinTable()
            bintable.extname("CIF")
            bintable.card("CIFVERSN", "1992a", "Version of CIF format")     
            cif.append(bintable)
            
            table = cif.table("CIF")
            
            # Copied from cta_root2caldb.py
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
            if caldir == self.cal_dir:
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
            table["TELESCOP"][row] = self.m_mission
            table["INSTRUME"][row] = self.m_instrument
            table["DETNAM"][row]   = "NONE"
            table["FILTER"][row]   = "NONE"
            table["CAL_DEV"][row]  = "ONLINE"
            table["CAL_CLAS"][row] = "BCF"
            table["CAL_DTYP"][row] = "DATA"
            table["CAL_VSD"][row]  = today
            table["CAL_VST"][row]  = time.split(".")[0]
            table["REF_TIME"][row] = 51544.0
            table["CAL_QUAL"][row] = 0
            table["CAL_CBD"][row]  = "NAME("+self.m_rspname+")"
            table["CAL_DATE"][row] = today.replace("-","/")[2:]
            table["CAL_DIR"][row]  = self.cal_dir
            table["CAL_FILE"][row] = self.m_outfile

        # Add effective area information
        row = row_index-4
        table["CAL_CNAM"][row]  = "EFF_AREA"
        table["CAL_DESC"][row]  = self.m_instrument+" effective area"

        # Add point spread function information
        row = row_index-3
        table["CAL_CNAM"][row]  = "RPSF"
        table["CAL_DESC"][row]  = self.m_instrument+" point spread function"
        
        # Add energy dispersion information
        row = row_index-2
        table["CAL_CNAM"][row]  = "EDISP"
        table["CAL_DESC"][row]  = self.m_instrument+" energy dispersion"

        # Add background information
        row = row_index-1
        table["CAL_CNAM"][row]  = "BGD"
        table["CAL_DESC"][row]  = self.m_instrument+" background"
        
        # Return CIF FITS file
        return cif
        
    def make_dirs(self):

        # Create gammalib Caldb        
        caldb = gammalib.GCaldb()
        
        # Set calibration directory 
        self.cal_dir  = "data"
        self.cal_dir += "/"+self.m_mission.lower()
        self.cal_dir += "/"+self.m_instrument.lower()
        self.cal_dir += "/bcf/"+self.m_rspname
        
        # Set absolute path
        self.base_dir = caldb.rootdir() +"/data"
        self.base_dir += "/"+self.m_mission.lower()
        self.base_dir += "/"+self.m_instrument.lower()
        
        # Set directory for irf file
        self.rsp_dir = caldb.rootdir() + "/" + self.cal_dir
        
        # Create directories and log information
        if not os.path.isdir(self.rsp_dir):
            if self.logExplicit():
                self.log(gammalib.parformat("Directory"))
                self.log(self.rsp_dir)
                self.log("\n")
            os.makedirs(self.rsp_dir)
        else:
            if self.logExplicit():
                self.log(gammalib.parformat("Directory (existing)"))
                self.log(self.rsp_dir)
                self.log("\n")

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
            self.log("\n\n")
            self.log.header1("Creating directory structure")
        
        # Create directory structure
        self.make_dirs()

        # Logging
        if self.logTerse():
            self.log("\n")
            self.log.header1("Updating calibration database")
         
        # Create/update calibration database   
        self.caldb_inx = self.make_caldb_index()
        
        # Logging
        if self.logTerse():
            self.log("\n")
            self.log.header1("Creating response file")
            
        # Create IRF file
        self.irf_fits = self.make_irf_file()
                
        # Return
        return

    def save(self):

        # Save caldb index file
        self.caldb_inx.save(self.m_clobber)

        # Save response file
        self.irf_fits.saveto(self.rsp_dir + "/" + self.m_outfile, self.m_clobber)
        
        # Return
        return
        
    def execute(self):
        """
        Execute the script.
        """
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
    """
    Generates caldb entry.
    """
    # Create instance of application
    app = csobs2caldb(sys.argv)
    
    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
