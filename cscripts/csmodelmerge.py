#!/usr/bin/env python
# ==========================================================================
# Merge model definition XML files
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
import glob


# ================== #
# csmodelmerge class #
# ================== #
class csmodelmerge(ctools.cscript):
    """
    Merge model definition XML files.

    An arbitrary number of model definition XML files will be merged into
    a single model definition XML file.    
    """

    # Constructors and destructors
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = "csmodelmerge"
        self._version = "1.1.0"

        # Initialise class members
        self._files      = None
        self._models     = gammalib.GModels()
        self._outmodel   = "NONE"

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
        Get parameters from parfile.
        """
        # Get input models string
        inmodels = self["inmodels"].string()
        
        # Handle ascii files
        if "@" == inmodels[0]:
            self._files = open(inmodels.replace("@","")).read().splitlines()  
            
        # Handle wild card strings
        elif "*" in inmodels:
            self._files = glob.glob(inmodels)
        
        # Handle space separated list
        elif " " in inmodels:
            self._files = inmodels.split(" ")
        
        # Handle semi-colon separated list
        elif ";" in inmodels:
            self._files = inmodels.split(";")
            
        # Throw exception if input models cannot be decoded
        else:
            msg = "Parameter \"inmodels\" must contain either an @ASCII "\
                  "file, a semi-colon-separated or whitespace-separated "\
                  "list of files or a wildcard string."
            raise RuntimeError(msg) 
        
        # Read ahead output filename
        if self._read_ahead():
            self.outmodel = self["outmodel"].filename()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
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
        
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Merge models")

        # Initialise model container
        self._models = gammalib.GModels()
        
        # Loop over model files
        for file in self._files:

            # Construct container from XML file
            models = gammalib.GModels(file)
            
            # Log number of models to add
            if self._logTerse():
                nmodels = models.size()
                if nmodels == 0:
                    self._log(gammalib.parformat("Add no model from file"))
                elif nmodels == 1:
                    self._log(gammalib.parformat("Add 1 model from file"))
                else:
                    self._log(gammalib.parformat("Add %d models from file" % 
                                                 nmodels))
                self._log(file)
                self._log("\n")
            
            # Extend model container by adding all models in the model file
            self._models.extend(models)
                
        # Log total number of models
        if self._logTerse():
            self._log(gammalib.parformat("Models after merging"))
            self._log(self._models.size())
            self._log("\n")

        # Return
        return

    def save(self):
        """ 
        Save model definition XML file
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save models")

        # Get output filename in case it was not read ahead
        self._outmodel = self["outmodel"].filename()
        
        # Log filename
        if self._logTerse():
            self._log(gammalib.parformat("Model definition XML file"))
            self._log(self._outmodel.url())
            self._log("\n")

        # Save models
        self._models.save(self._outmodel)
        
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

        # Save ds9 file if required
        self.save()

        # Return
        return    
        

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csmodelmerge(sys.argv)
    
    # Execute application
    app.execute()
