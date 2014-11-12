#! /usr/bin/env python
# ==========================================================================
# Script for spectral data points generation.
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

import os


# ============== #
# cstsdist class #
# ============== #
class csspec(ctools.ctool):
    """
    This class implements the creation of spectral points. It derives
    from the ctools::ctool class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool. 
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        self.name    = "csspec"
        self.version = "0.1.0"
        #super(ctools.ctool).__init__(self.name,self.version,argv)
        # Set name
#         self.name    = "csspectrum"
#         self.version = "0.1.0"
        
        # Initialise some members
        self.obs       = None 
              
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
            pars.append(gammalib.GApplicationPar("srcname","s","a","Crab","","","Source name"))
            pars.append(gammalib.GApplicationPar("outfile","f","h","spectrum.fits","","","Output file name"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","$GAMMALIB/share/caldb/cta","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","cta_dummy_irf","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("emin","r","a","0.1","0.0","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","a","100.0","0.0","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","20","","","Number of energy bins"))
            pars.append(gammalib.GApplicationPar("ebinalg","s","h","LOG","LIN|LOG|FILE","","Binning algorithm"))
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
            self.models = gammalib.GModels(self["srcmdl"].filename())
            
        self.m_srcname   = self["srcname"].string()
        self.m_outfile   = self["outfile"].filename()


        # not sure how to call get_ebounds as private member of ctool        
        #self.ebounds = ctools.ctool.get_ebounds()
        self.m_emin      = self["emin"].real()
        self.m_emax      = self["emax"].real()
        self.m_enumbins  = self["enumbins"].integer()
        #self.m_ebinalg   = self["ebinalg"].string()
        
        self.ebounds = gammalib.GEbounds(self.m_enumbins,gammalib.GEnergy(self.m_emin,"TeV"),gammalib.GEnergy(self.m_emax,"TeV"),True)

             
        # Set some fixed parameters
        self.m_log   = False # Logging in client tools
        self.m_debug = False # Debugging in client tools
        self.m_clobber = self["clobber"].boolean()
        
        self.parname =self.models[self.m_srcname].spectral()[0].name()
  
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
        self.save(self.m_outfile, self.m_clobber)
        
        # Return
        return

    def analyse_bin(self, emin,emax):
        #print "1",self.models.size()
        select = ctools.ctselect(self.obs)  
        select["emin"].real(emin)  
        select["emax"].real(emax) 
        select["ra"].real(-1.0)
        select["dec"].real(-1.0)
        select["rad"].real(-1.0)
        select["tmin"].real(0.0)
        select["tmax"].real(0.0)
        select.run()
        #print "2",self.models.size()
        like = ctools.ctlike(select.obs())
        like.run()
        optmodels = like.obs().models()
        #print "3",self.models.size()
        return optmodels.copy()

    def run(self):
        """
        Run the script.
        """
        
        # Switch screen logging on in debug mode
        if self.logDebug():
            self.log.cout(True)

        # Get parameters
        self.get_parameters()

        
        for model in self.models:
            for par in model:
                par.fix()
                  
        self.models[self.m_srcname][self.parname].free()
        self.models[self.m_srcname].tscalc(True)
        self.obs.models(self.models)
        
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
        
        binned=False
        
        
        for i in range(self.ebounds.size()):
            
            emin = self.ebounds.emin(i).TeV()
            emax = self.ebounds.emax(i).TeV()
            emean = self.ebounds.elogmean(i)
            optmodels = self.analyse_bin(emin,emax)
            #print optmodels
            #if optmodels.size()>0:
            optmodel = optmodels[self.m_srcname]
            ts = optmodel.ts()
            
            flux = optmodel.spectral().eval(emean,gammalib.GTime())
            rel_error = optmodel[self.parname].error()/optmodel[self.parname].value()
            
            flux_error = flux * rel_error
            
            print emean,flux,flux_error,ts
        
            if binned:
                print "not implemented"

        self.obs.models(self.models)
        
        # Return
        return
    
    def save(self,outfile,clobber):
        return

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates residual count map.
    """
    # Create instance of application
    app = csspec(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
    