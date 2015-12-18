#!/usr/bin/env python
# ==========================================================================
# Merges sliced TS map files into one map
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
import glob
import os


# ================== #
# cstsmapmerge class #
# ================== #
class cstsmapmerge(ctools.cscript):
    """
    This class merges TS map files into one TS map
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "cstsmapmerge"
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
            pars.append(gammalib.GApplicationPar("inmaps","s","a","tsmap.fits","","","Input TS map FITS files"))
            pars.append(gammalib.GApplicationPar("outmap","f","a","NONE","","","Output TS map FITS file"))
            pars.append(gammalib.GApplicationPar("overwrite","b","h","yes","","","Overwrite previously filled bins?"))
            pars.append(gammalib.GApplicationPar("delete","b","h","no","","","Delete merged files?"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","cstsmapmerge.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get input models
        self.inmodels = self["inmaps"].string()
        
        # Interpet inmodels string
        # Handle ascii files
        if "@" == self.inmodels[0]:
            self.files = open(self.inmodels.replace("@","")).read().splitlines()  
            
        # Handle wild card strings
        elif "*" in self.inmodels:
            self.files = glob.glob(self.inmodels)
        
        # Handle space-separated list
        elif " " in self.inmodels:
            self.files = self.inmodels.split(" ")
        
        # Handle comma-separated list
        elif "," in self.inmodels:
            self.files = self.inmodels.split(",")
            
        # Throw exception if input models cannot be decoded
        else:
            msg = "Parameter \"inmmaps\" must contain either an @ASCII file,"+\
                  " a comma-separated or whitespace-separated list of files"+\
                  " or a wildcard string"
            raise RuntimeError(msg) 

        # Check number of files. We need at least two files.
        if len(self.files) <= 1:
            msg = "Need at least two files to start merging, "+\
                  str(len(self.files))+" files given."
            raise RuntimeError(msg) 
        
        # Get output file
        self.outmap = self["outmap"].filename()
        
        # Get other parameters
        self.overwrite = self["overwrite"].boolean()
        self.delete = self["delete"].boolean()
        
        # Return
        return
 
    def init_map(self, fitsfile):
        """
        """
        # Set filename
        self.in_filename = fitsfile
        
        # Open FITS file
        fits           = gammalib.GFits(fitsfile)
        self.tsmap     = gammalib.GSkyMap()
        self.tsmap.read(fits[0])
        self.statusmap = gammalib.GSkyMap()
        self.statusmap.read(fits["STATUS MAP"])
        
        # Get other maps 
        self.maps     = []
        self.mapnames = []
        
        # Loop over extensions
        for hdu in fits:
            
            # Leave out primary and status extension
            if hdu.extname() != "IMAGE" and hdu.extname() != "STATUS MAP":
                
                # Add present maps
                skymap = gammalib.GSkyMap()
                skymap.read(hdu)
                self.maps.append(skymap)
                self.mapnames.append(hdu.extname())
                
        # Close FITS file       
        fits.close()
        
        # Return
        return

    # Add map instances to each other. In this way the TS maps get merged  
    def add(self, fitsfile):
        
        # open FITS file
        fits          = gammalib.GFits(fitsfile)
        add_tsmap     = gammalib.GSkyMap()
        add_tsmap.read(fits[0])
        add_statusmap = gammalib.GSkyMap()
        add_statusmap.read(fits["STATUS MAP"])
        
        # Get other maps 
        add_maps = []
        
        # Loop over extensions
        for hdu in fits:
            
            # Leave out primary and status extension
            if hdu.extname() != "IMAGE" and hdu.extname() != "STATUS MAP":
                
                # Add present maps
                skymap = gammalib.GSkyMap()
                skymap.read(hdu)
                add_maps.append(skymap)
                
        # Close FITS file
        fits.close()
        
        # Compare size of other maps
        if not len(add_maps) == len(self.maps):
            self.error("Maps cannot be added: Maps from \""+self.in_filename+"\" do not match from \""+fitsfile+"\"")

        # Loop over bins    
        for i in range(self.tsmap.npix()):
            if add_statusmap[i] > 0.5:
                
                # Throw exception if this bin has already been computed
                if self.statusmap[i] > 0.5 and not self.overwrite:
                    msg = "Attempt to merge bin which apparently has already"+\
                          " been merged. File \""+fitsfile+"\" contains"+\
                          " already merged bins. Set hidden parameter"+\
                          " overwrite=yes to avoid this error."
                    raise RuntimeError(msg)
                
                # Copy TS values
                self.tsmap[i] = add_tsmap[i]
                
                # Copy status 
                self.statusmap[i] = add_statusmap[i]
                
                # Loop over maps and copy entries
                for j in range(len(self.maps)):
                    self.maps[j][i] = add_maps[j][i]
        
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
            
        if self.logTerse():
            self.log("\n")
            self.log.header1("Merging sliced TS map files")
            self.log("\n")
        
        # Initialise file to start with
        file0 = ""
        
        self.added_files = []
        
        # Test files for the entry status map
        # use the first one to appear useful
        for fitsfile in self.files:    
            if not gammalib.is_fits(fitsfile):
                self.log("Skipping file \""+fitsfile+" (not a FITS file)\n")
            fits = gammalib.GFits(fitsfile)
            
            if fits.contains("STATUS MAP"):
                self.log("Using file \""+fitsfile+"\" as initialising map\n")
                file0 = fitsfile
                break
            else:
                if self.logExplicit():
                    self.log("File \""+fitsfile+"\" has not \"STATUS MAP\" extension\n")
        
        # Print an error message if no suitable file was provided
        if file0 == "":
            if self.logTerse():
                self.log("None of the provided file seems to be a sliced ts map file (no \"STATUS MAP\" extension)\n") 
        
        # Otherwise merge files
        else:
            
            # copy file list
            workfiles = self.files
            
            # Remove entry which will be used to initalise the map
            workfiles.remove(file0)
            
            # Initialise map from first file
            self.init_map(file0) 
            
            # Append to added files 
            self.added_files.append(file0)
            
            # Loop over files
            for fitsfile in workfiles:
                
                # Skip if file is not FITS
                if not gammalib.is_fits(fitsfile):
                    if self.logTerse():
                        self.log("Skipping file \""+fitsfile+" (not a FITS file)\n") 
                    continue   
                
                # Skip if file does not contain status map
                fits = gammalib.GFits(fitsfile)
                if not fits.contains("STATUS MAP"):
                    if self.logTerse():
                        self.log("Skipping file \""+fitsfile+" (not a FITS file)\n")    
                    continue
                
                # Logging
                if self.logTerse():
                    self.log("Adding TS map slice \""+fitsfile+"\"\n")  
                    
                # Add file 
                self.add(fitsfile)
                
                # Append to added files
                self.added_files.append(fitsfile)
                
        # Return
        return

    def save(self):
        """ 
        Save merged tsmap file and remove slices if requested
        """
        
        # Logging
        if self.logExplicit():
            self.log("\n")
            self.log.header1("Saving merged TS map")
            self.log("\n") 
        
        # Create FITS file
        fits = gammalib.GFits()
        
        # Write TS map into primary
        self.tsmap.write(fits)

        # Loop over maps and write them to fits
        for i in range(len(self.maps)):
            self.maps[i].write(fits)
        
        # Set map names as extensions
        for i in range(len(self.mapnames)):   
            fits[i+1].extname(self.mapnames[i])

        # Check if map is fully done
        done = True
        for pix in self.statusmap:
            if pix < 0.5:
                done = False
                break
        
        # Write Status map if we are not done yet
        if not done:
            self.statusmap.write(fits)
            fits[fits.size()-1].extname("STATUS MAP")
        
        # Save FITS file
        fits.saveto(self.outmap, self.clobber())
        if self.logExplicit():
            self.log("Saved output map to \""+self.outmap+"\"\n")
        if self.delete:
            if self.logTerse():
                self.log("\n")
                self.log.header1("Deleting slices")
                self.log("\n")
            for filename in self.added_files:
                if self.logExplicit():
                    self.log("Deleting file \""+filename+"\"\n")
                os.remove(filename)
        
        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()

        # Save ds9 file if required
        self.save()

        # Return
        return    
        

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Merges ts map slices.
    """
    # Create instance of application
    app = cstsmapmerge(sys.argv)
    
    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
