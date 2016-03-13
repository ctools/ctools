#!/usr/bin/env python
# ==========================================================================
# Find observations from an IACT data store
#
# Copyright (C) 2016 Michael Mayer
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


# =============== #
# csfindobs class #
# =============== #
class csfindobs(ctools.cscript):
    """
    Find observations from an IACT data store.
    """

    # Constructors and destructors
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name     = "csfindobs"
        self._version  = "1.1.0"
        self._verbose  = False
        self._datapath = os.getenv("VHEFITS","")
        
        # Make sure that parfile exists
        self._parfile()

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self._name, self._version)
        elif len(argv) ==1:
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
        #  Write separator into logger
        if self._logTerse():
            self._log("\n")
        
        # Return
        return

    def _parfile(self):
        """
        Check if parfile exists. If parfile does not exist then create a
        default parfile. This kluge avoids shipping the cscript with a
        parfile.
        """
        # Set parfile name
        parfile = self._name+".par"
        
        try:
            pars = gammalib.GApplicationPars(parfile)
        except:
            # Signal if parfile was not found
            sys.stdout.write("Parfile "+parfile+" not found. Create default parfile.\n")
            
            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("datapath","s","a",self._datapath,"","","Path were data are located"))       
            pars.append(gammalib.GApplicationPar("prodname","s","a","prod-name","","","Name of FITS production (Run csiactdata to view your options)"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","runlist.lis","","","Runlist outfile"))
            pars.append(gammalib.GApplicationPar("ra","r","a","83.6331","","","Right ascension"))
            pars.append(gammalib.GApplicationPar("dec","r","a","22.01","","","Declination"))
            pars.append(gammalib.GApplicationPar("rad","r","a","2.5","","","Search radius"))
            pars.append(gammalib.GApplicationPar("min_qual","i","h","0","0|1|2","","Minimum data quality (0=perfect, 1=ok, 2=bad)"))
            pars.append(gammalib.GApplicationPar("expression","s","h","NONE","","","Additional expression"))
            pars.append(gammalib.GApplicationPar("master_indx","s","h","master.json","","","Name of master index file"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csfindobs.log","","","Log filename"))
            pars.save(parfile)
        
        # Return
        return
        
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters
        if self._datapath == "":
            self._datapath = self["datapath"].string()
        
        # Expand environment
        self._datapath = gammalib.expand_env(self._datapath)
        
        # Get production name
        self._prodname = self["prodname"].string()
        
        # Master index file name
        self._master_indx = self["master_indx"].string()
        
        # Intiialise flag if spatial selection is required
        self._select_radec = True
        
        # Initialise invalid radius
        self._radius = 0.0
        
        # Check for vailidity of spatial parameters
        if (self["ra"].is_valid() and
            self["dec"].is_valid() and 
            self["rad"].is_valid()):
            
            # Read spatial parameters
            self._ra     = self["ra"].real()
            self._dec    = self["dec"].real()
            self._radius = self["rad"].real()
            
        else:
            self._select_radec = False
        
        # Check Radius for validity
        if self._radius <= 0.0:
            self._select_radec = False
            
        # Read other parameters
        self._quality    = self["min_qual"].integer()
        self._expression = self["expression"].string()
        
        # Check for validity of expression
        if self._expression == "NONE" or self._expression == "INDEF":
            self._expression = ""
        
        # Read outfile
        self._outfile = self["outfile"].filename()
        
        # Open master index file and look for prodname 
        master_file = os.path.join(self._datapath, self._master_indx)
        if not os.path.isfile(master_file):
            raise RuntimeError("FITS data store not available. "+
                               "No master index file found at \""+
                               master_file+"\". Make sure the file "+
                               "is copied from the server and your "+
                               "datapath is set correctly.")

        # Open and load JSON file
        json_data = open(master_file).read()
        data      = json.loads(json_data)    
        if not "datasets" in data:
            raise RuntimeError("Key \"datasets\" not available in "+
                               "master index file.")

        # Get configurations
        configs = data["datasets"]

        # Initialise obs index file
        self._obs_index = ""

        # Get HDUs
        for config in configs:
            if self._prodname == config["name"]:
                self._obs_index = str(os.path.join(self._datapath,
                                                     config["obsindx"]))
                break

        # Check HDUs
        if self._obs_index == "":
            raise RuntimeError("FITS data store \""+self._prodname+
                               "\" not available. Run csiactdata to get "+
                               "a list of available storage names")
        filename = gammalib.GFilename(self._obs_index+"[OBS_INDEX]")
        if not filename.is_fits():
            raise RuntimeError("Observation index file \""+
                               self._obs_index+
                               "[OBS_INDEX]\" for FITS data store \""+
                               self._prodname+
                               "\" not available. "+
                               "Check your master index file or run "+
                               "csiactdata to get a list of available "+
                               "storage names.")

        # Set other members
        self._debug   = False # Debugging in client tools
        self._clobber = self["clobber"].boolean()
  
        # Return
        return

    def obs_ids(self):
        """ 
        Return OBS IDs
        """
             
        # Return
        return self.runs
    
    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()
        
        # Save residual map
        self.save(self._outfile)
        
        # Return
        return

    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()
        
        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")
        
        # Initialise run list
        self._runs = []
        
        # Initialise selction expression
        expr = ""
        
        # Add spatial expression if possible
        if self._select_radec:
            expr += "ANGSEP("+str(self._ra)+","+str(self._dec)+ \
                    ",RA_PNT,DEC_PNT)<="+str(self._radius)
        
        # Add connector if expression is empty
        if len(expr):
            expr += "&&"
        
        # Add quality expression
        expr += "QUALITY<="+str(self._quality)
        
        # Add user expression      
        if self._expression == "NONE" or len(self._expression) > 0:
            expr += "&&"+self._expression   
            
        # Initialise filename   
        openfile = self._obs_index+"[OBS_INDEX]["+expr+"]"    
        
        # Open file
        fits = gammalib.GFits(openfile)

        # Check for extistence of rows
        if fits["OBS_INDEX"].nrows() > 0:
            
            # Get HDU
            obs_indx = fits["OBS_INDEX"]
            
            # Read observation IDs into array
            for i in range(obs_indx.nrows()):
                self._runs.append(obs_indx["OBS_ID"][i])

        # Dump expression
        if self._logNormal():
            self._log("\n")
            self._log.header3("Expression")
            self._log(expr)
            self._log("\n")

        # Dump number of observations
        if self._logTerse():
            self._log("\n")
            self._log.header3("Observations")
            self._log("Found "+str(len(self._runs))+" observations")
            self._log("\n")
        
        # Close file
        fits.close()

        # Return
        return
    
    def save(self, outfile):
        """
        Save.
        """
        # Check for clobber
        if outfile.exists() and not self._clobber():
            self._log("File \""+outfile+"\" already exists. ")
            self._log("Set clobber=yes to overwrite file.")
            self._log("\n")
        else:
            # Write runs to file
            f = open(repr(outfile),"w")
            for run in self._runs:
                f.write(str(run)+" \n")
            f.close()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csfindobs(sys.argv)
    
    # Open logfile
    app._logFileOpen()
    
    # Execute application
    app.execute()
