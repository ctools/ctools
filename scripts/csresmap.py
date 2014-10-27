#! /usr/bin/env python
# ==========================================================================
# Residual map generation script.
#
# Copyright (C) 2014 Michael Mayer
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
import ctools
import gammalib
import sys
import tempfile
import os
import math


# ============== #
# cstsdist class #
# ============== #
class csresmap(gammalib.GApplication):
    """
    This class implements the creation of a residual map. It derives
    from the GammaLib::GApplication class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool. 
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "csresmap"
        self.version = "0.1.0"
        
        # Initialise some members
        self.obs       = None 
        self.algorithm = "SUB"
        self.resmap    = None
              
        # Initialise some members
        if isinstance(argv[0],gammalib.GObservations):
            self.obs = argv[0]
            argv = argv[1:]
        else:      
            self.obs      = gammalib.GObservations()
            self.obs.clear()   
        self.outfile = ""
        
        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            gammalib.GApplication.__init__(self, self.name, self.version)
        elif len(argv) ==1:
            gammalib.GApplication.__init__(self, self.name, self.version, *argv)
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
        #  Write separator into logger
        if self.logTerse():
            self.log("\n")
        
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
            # Signal if parfile was not found
            sys.stdout.write("Parfile "+parfile+" not found. Create default parfile.\n")
            
            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("infile","f","a","obs.xml","","","Observation definition"))
            pars.append(gammalib.GApplicationPar("srcmdl","f","a","crab.xml","","","model file name"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","resmap.fits","","","Output file name"))
            pars.append(gammalib.GApplicationPar("binfile","f","h","binned.fits","","","Output binned file name"))
            pars.append(gammalib.GApplicationPar("modfile","f","h","model.fits","","","Output model file name"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","$GAMMALIB/share/caldb/cta","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","cta_dummy_irf","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("emin","r","h","0.1","0.0","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","h","100.0","0.0","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","h","20","","","Number of energy bins"))
            pars.append(gammalib.GApplicationPar("ebinalg","s","h","LOG","LIN|LOG|FILE","","Binning algorithm"))
            pars.append(gammalib.GApplicationPar("coordsys","s","a","CEL","CEL|GAL","","Coordinate System"))
            pars.append(gammalib.GApplicationPar("proj","s","a","TAN","AIT|AZP|CAR|MER|STG|TAN","","Projection method e.g. AIT|AZP|CAR|MER|STG|TAN"))
            pars.append(gammalib.GApplicationPar("xref","r","a","83.6331","","","First coordinate of image center in degrees (RA or galactic l)"))
            pars.append(gammalib.GApplicationPar("yref","r","a","22.01","","","Second coordinate of image center in degrees (DEC or galactic b)"))
            pars.append(gammalib.GApplicationPar("nxpix","i","h","200","","","Size of the X axis in pixels"))
            pars.append(gammalib.GApplicationPar("nypix","i","h","200","","","Size of the Y axis in pixels"))
            pars.append(gammalib.GApplicationPar("binsz","r","h","0.05","","","Pixel size (deg/pixel)"))
            pars.append(gammalib.GApplicationPar("algorithm","s","h","SUBDIV","SUB|SUBDIV|SUBDIVSQRT","","Residualmap algorithm"))
            pars.append_standard()
            pars.save(parfile)
        
        # Return
        return
        
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters
        
        # Set observation if not done before
        if self.obs.size() == 0:
            
            obsfile = self["infile"].filename()
            try: 
                
                self.obs = gammalib.GObservations(obsfile)
            except:
                self.obs.clear()
                observation = gammalib.GCTAObservation()
                observation.load(obsfile)
                
                self.m_irf   = self["irf"].string()
                self.m_caldb = self["caldb"].string()
                
                caldb = gammalib.GCaldb()
                if os.path.isdir(self.m_caldb):
                    caldb.rootdir(self.m_caldb)    
                else:
                    caldb.open("cta",self. m_caldb)

                observation.response(self.m_irf, caldb);
                self.obs.append(observation)
            
        # Check for models in the container
        # read models from file if there were none
        # Note that the models should be optimised for 
        # meaningful residual map    
        if self.obs.models().size() == 0:
            self.obs.models(self["srcmdl"].filename())
        self.m_outfile   = self["outfile"].filename()
        self.m_modfile   = self["modfile"].filename()
        self.m_binfile   = self["binfile"].filename()
        self.m_xref      = self["xref"].real()
        self.m_yref      = self["yref"].real()
        self.m_emin      = self["emin"].real()
        self.m_emax      = self["emax"].real()
        self.m_enumbins  = self["enumbins"].integer()
        self.m_ebinalg   = self["ebinalg"].string()
        self.m_coordsys  = self["coordsys"].string()
        self.m_proj      = self["proj"].string()
        self.m_nxpix     = self["nxpix"].integer()
        self.m_nypix     = self["nypix"].integer()
        self.m_binsz     = self["binsz"].real()
        self.m_algorithm = self["algorithm"].string()
             
        # Set some fixed parameters
        self.m_log   = False # Logging in client tools
        self.m_debug = False # Debugging in client tools
        self.m_clobber = self["clobber"].boolean()
  
        # Return
        return
    
    def models(self, models):
        """
        Set model.
        """
        # Copy models
        self.obs.models(models.clone())
    
        # Return
        return
        
    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()

        # Save residual map
        self.resmap.save(self.m_outfile, self.m_clobber)
        
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
        
        #  Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")
        
        # Write observation into logger
        if self.logTerse():
            self.log("\n")
            self.log.header1("Observation")
            self.log(str(self.obs))
            self.log("\n")

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Generate binned map")

        # Create countsmap
        bin = ctools.ctbin(self.obs)
        bin["debug"].boolean(True)
        bin["nxpix"].integer(self.m_nxpix)
        bin["nypix"].integer(self.m_nypix)
        bin["proj"].string(self.m_proj)
        bin["coordsys"].string(self.m_coordsys)
        bin["xref"].real(self.m_xref)
        bin["yref"].real(self.m_yref)
        bin["enumbins"].integer(self.m_enumbins)
        bin["ebinalg"].string(self.m_ebinalg)
        bin["emin"].real(self.m_emin)
        bin["emax"].real(self.m_emax)
        bin["binsz"].real(self.m_binsz)
        bin["clobber"].boolean(self.m_clobber)
        bin.run()

        # Store counts map as residual map. Note that we need a
        # special construct here to avoid memory leaks. This seems
        # to be a SWIG feature as SWIG creates a new object when
        # calling bin.cube()
        #residualmap = bin.cube().map()
        counts_cube = bin.cube()
        self.resmap = counts_cube.map().copy()
        self.resmap.stack_maps()
 
        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Generate model map")
    
        # Create model map
        model = ctools.ctmodel(self.obs)
        model["debug"].boolean(True)
        model.cube(bin.cube())
        model["clobber"].boolean(self.m_clobber)
        model.run()

        # Get model map into GSkymap object
        modelmap = model.cube().map()
        modelmap.stack_maps()
        
        # Continue calculations depending on given algorithm
        if self.m_algorithm == "SUB":
            
            # Subtract maps 
            self.resmap -= modelmap
        
        elif self.m_algorithm == "SUBDIV":
            
            # Subtract and divide by model map
            self.resmap -= modelmap
            #self.resmap /= modelmap   # Python 3.x does not like this !!!
            for pixel in modelmap:
                if pixel != 0.0:
                    pixel = 1.0/pixel
            self.resmap *= modelmap
            
        elif self.m_algorithm == "SUBDIVSQRT":

            # subtract and divide by sqrt of model map
            self.resmap -= modelmap
            #self.resmap /= modelmap.sqrt()   # Python 3.x does not like this !!!
            for pixel in modelmap:
                if pixel != 0.0:
                    pixel = 1.0/math.sqrt(pixel)
            self.resmap *= modelmap
            
        else:
            
            # Raise error if algorithm is unkown
            raise TypeError("Algorithm \""+self.m_algorithm+"\" not known")
        
        # Return
        return
    

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates residual count map.
    """
    # Create instance of application
    app = csresmap(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
    