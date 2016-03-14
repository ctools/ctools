#!/usr/bin/env python
# ==========================================================================
# Inspection of IACT data storage
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
import json


# ================ #
# csiactdata class #
# ================ #
class csiactdata(ctools.cscript):
    """
    This script inspects the available FITS data storage on the user machine
    and writes information about the available FITS configurations into a log
    file or prints it on the screen. These information can be used as input
    for 'csiactobs'.
    """

    # Constructors and destructors
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name      = "csiactdata"
        self._version   = "1.1.0"
        self._datapath  = os.getenv("VHEFITS","")
        self._prodnames = []

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
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        
        # Get datapath if not already set
        if self._datapath == "":
            self._datapath = self["datapath"].string()
        
        # Expand environment
        self._datapath = gammalib.expand_env(self._datapath)
        
        # Get filename of master index file
        self._master_indx = self["master_indx"].string()
        self._master_file = os.path.join(self._datapath, self._master_indx)
        
        # Check for presence of master index file
        if not os.path.isfile(self._master_file):
            raise RuntimeError("Master index file \""+self._master_file+"\" not found. Use hidden parameter \"master_indx\" to specifiy a different filename.")
        
        # Return
        return

    # Public methods
    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on
        self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")
        
        self._log("\n")
        self._log.header1("Data storage entries")
        self._log("\n")
        self._log.parformat("Master index file")
        self._log(self._master_file)
        self._log("\n")
        self._log("\n")
        
        # Open and load JSON file
        json_data = open(self._master_file).read()
        data      = json.loads(json_data)    
        configs   = data["datasets"]
        
        # Initialise array for available names
        self._prodnames = []
        
        # Loop over configs and display unavailable storage first        
        for config in configs:
            
            # Create hdu and obs index files
            hdu = os.path.join(self._datapath, config["hduindx"])
            obs = os.path.join(self._datapath, config["obsindx"])
            
            filename_hdu = gammalib.GFilename(str(hdu)+"[HDU_INDEX]")
            filename_obs = gammalib.GFilename(str(obs)+"[OBS_INDEX]")
            
            # Check if index files are available
            if not (filename_hdu.is_fits() and filename_obs.is_fits()):
                self._log.parformat(str(config["name"]))
                self._log("Not available")
                self._log("\n")
        
        
        self._log("\n")
        self._log.header2("Available data configs")
        
        # Loop over configs and log available configs       
        for config in configs:       
            
            # Create hdu and obs index files
            hdu = os.path.join(self._datapath, config["hduindx"])
            obs = os.path.join(self._datapath, config["obsindx"])
            
            filename_hdu = gammalib.GFilename(str(hdu)+"[HDU_INDEX]")
            filename_obs = gammalib.GFilename(str(obs)+"[OBS_INDEX]")
            
            # Check if index files are available
            if (filename_hdu.is_fits() and filename_obs.is_fits()):
                
                # Log data information  
                self._log("\n")
                self._log.header3(str(config["name"]))
                
                # Print important informations first
                self._log.parformat("Name")
                self._log(str(config["name"]))
                self._log("\n")  
                self._log.parformat("Observation index")
                self._log(str(config["obsindx"]))
                self._log("\n") 
                self._log.parformat("HDU index")
                self._log(str(config["hduindx"]))
                self._log("\n") 
                
                # Append to available names
                self._prodnames.append(str(config["name"]))
                
                # Print additional information
                for key in config:
                    if key in ["name","hduindx","obsindx"]:
                        continue
                    self._log.parformat(str(key))
                    self._log(str(config[key]))
                    self._log("\n")  

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

        # Return
        return


    def names(self):
        """ 
        Return available FITS production names
        """
        
        # Return 
        return self._prodnames


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csiactdata(sys.argv)

    # Execute application
    app.execute()
