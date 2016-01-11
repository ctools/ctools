#!/usr/bin/env python
# ==========================================================================
# Inspection of IACT data storage
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
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name     = "csiactdata"
        self.version  = "1.0.0"
        self.datapath = os.getenv("VHEFITS","")

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
            # Signal if parfile was not found
            sys.stdout.write("Parfile "+parfile+" not found. Create default parfile.\n")

            # Create default parfile
            pars = gammalib.GApplicationPars()    
            pars.append(gammalib.GApplicationPar("datapath","s","a",self.datapath,"","","Path were data is located"))     
            pars.append(gammalib.GApplicationPar("master_indx","s","h","master.json","","","Name of master index file"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csiactdata.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        
        # Get datapath if not already set
        if self.datapath == "":
            self.datapath = self["datapath"].string()
        
        # Expand environment
        self.datapath = gammalib.expand_env(self.datapath)
        
        # Get filename of master index file
        self.m_master_indx = self["master_indx"].string()
        self.m_master_file = os.path.join(self.datapath, self.m_master_indx)
        
        # Check for presence of master index file
        if not os.path.isfile(self.m_master_file):
            raise RuntimeError("Master index file \""+self.m_master_file+"\" not found. Use hidden parameter \"master_indx\" to specifiy a different filename.")
        
        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()

        # Return
        return

    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on
        self.log.cout(True)

        # Get parameters
        self.get_parameters()

        # Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")
        
        self.log("\n")
        self.log.header1("Data storage entries")
        self.log("\n")
        self.log.parformat("Master index file")
        self.log(self.m_master_file)
        self.log("\n")
        self.log("\n")
        
        # Open and load JSON file
        json_data = open(self.m_master_file).read()
        data      = json.loads(json_data)    
        configs   = data["datasets"]
        
        # Loop over configs and display unavailable storage first        
        for config in configs:
            
            # Create hdu and obs index files
            hdu = os.path.join(self.datapath, config["hduindx"])
            obs = os.path.join(self.datapath, config["obsindx"])
            
            # Check if index files are available
            if not (gammalib.is_fits(str(hdu)) and gammalib.is_fits(str(obs))):
                self.log.parformat(str(config["name"]))
                self.log("Not available")
                self.log("\n")
        
        # Loop over configs and log available configs       
        for config in configs:       
            
            # Create hdu and obs index files
            hdu = os.path.join(self.datapath, config["hduindx"])
            obs = os.path.join(self.datapath, config["obsindx"])
            
            # Check if index files are available
            if gammalib.is_fits(str(hdu)) and gammalib.is_fits(str(obs)):   
                
                # Log data information  
                self.log("\n")
                self.log.header3(str(config["name"]))
                for key in config:
                    self.log.parformat(str(key))
                    self.log(str(config[key]))
                    self.log("\n")  

        # Return
        return       

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Shows data storage information
    """
    # Create instance of application
    app = csiactdata(sys.argv)

    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
