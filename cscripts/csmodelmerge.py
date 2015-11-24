#!/usr/bin/env python
# ==========================================================================
# Merges model XML files into one file
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


# =============== #
# cshessobs class #
# =============== #
class csmodelmerge(ctools.cscript):
    """
    This class merges an arbitray number of model container into one output model XML file
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "csmodelmerge"
        self.version = "0.1.0"

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
        default parfile. This kluge avoids shipping the cscript with a parfile.
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
            pars.append(gammalib.GApplicationPar("inmodels","s","a","model.xml","","","Input model XML files"))
            pars.append(gammalib.GApplicationPar("outmodel","f","a","NONE","","","Output model file"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csmodelmerge.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get input models
        self.inmodels = self["inmodels"].string()
        
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
            msg = "Parameter \"inmodels\" must contain either an @ASCII file, a comma-separated or space-separated list of files or a wildcard string"
            raise gammalib.GException.invalid_argument("csmodelmerge::get_parameters", msg) 
        
        # Get output file
        self.outmodel = self["outmodel"].filename()

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
            self.log.header1("Merging Models")
            self.log("\n")
            
            
        # Initialise output model container
        self.models = gammalib.GModels()
        
        # Loop over model files
        for modelfile in self.files:
            
            # Logging
            self.log("Adding "+modelfile)
            self.log("\n")
            
            # Open models
            model_add = gammalib.GModels(modelfile)
            
            # Append models to container
            for model in model_add:
                self.models.append(model)
                
        # Return
        return

    def save(self):
        """ 
        Save pointings to ds9 region file if required
        """
        
        # Logging
        if self.logExplicit():
            self.log("\n")
            self.log.header1("Saving model file")
            self.log("\n") 
        
        # Save models
        self.models.save(self.outmodel)
        
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
    Merges model XML files
    """
    # Create instance of application
    app = csmodelmerge(sys.argv)
    
    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
