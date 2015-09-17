#!/usr/bin/env python
# ==========================================================================
# Generation of an H.E.S.S. observation definition file.
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
import glob
import math


# =============== #
# cshessobs class #
# =============== #
class cshessobs(ctools.cscript):
    """
    This class implements the creation of a observation xml file for HESS
    data analysis. This class is dedicated for use inside the HESS Collaboration,
    i.e. it can only be used if you have access to HESS data in FITS format.
    The HESSFITS environment has to be set correctly to use this script
    Please see the internal documentation for this purpose:
    https://hess-confluence.desy.de/confluence/display/HESS/HESS+Open+Source+Tools+-+HESS+data+and+IRFs+in+FITS+format
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "cshessobs"
        self.version = "0.2.0"
        self.datapath = ""
        self.m_erange = gammalib.GEbounds()

        # Check on existence of HESSFITS environment variable
        try:
            self.datapath = os.environ["HESSFITS"] 
        except KeyError:
            raise RuntimeError("Environment variable \"HESSFITS\" not set. " \
                               "Please set the environment variable \"HESSFITS\" to use this script")

        # Initialise some members
        self.inmodels = None
        self.outobs   = "obs.xml"

        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self.name, self.version)
        elif len(argv) ==1:
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

            # Check for folders and possible/available configurations
            chains = os.listdir(self.datapath)

            # Initialise folder dictionary
            folders = {}

            # Loop over possible analysis chains
            for chain in chains:

                # Omit hidden files/folders    
                if chain[0]=="." or not os.path.isdir(self.datapath+"/"+chain):
                    continue  

                # Initialise dictionary per chain
                folders[chain] = {}

                # Read available DST versions
                dstversions = os.listdir(self.datapath+"/"+chain)

                # Loop over DST versions
                for dstver in dstversions:

                    # Omit hidden files/folders    
                    if dstver[0] == "." or not os.path.isdir(self.datapath+"/"+chain+"/"+dstver):
                        continue

                    # Initialise available cut configurations
                    folders[chain][dstver] = []
                    configs = os.listdir(self.datapath+"/"+chain+"/"+dstver)

                    # Loop over configurations
                    for config in configs:   

                        # Omit hidden files/folders                   
                        if config[0] == "." or not os.path.isdir(self.datapath+"/"+chain+"/"+dstver+"/"+config):
                            continue

                        # Append configuration to corresponding dictionary
                        folders[chain][dstver].append(config)

            # Check for appropriate folder structure
            checked_folders = self.check_structure(folders)
            if not len(checked_folders):
                raise KeyError("Appropriate folder structure not found in $HESSFITS. The correct structure can be found on HESS internal documentation")

            # Create default parfile
            pars = gammalib.GApplicationPars()

            # Create parameters for each available config            
            avchains = ""
            for chain in checked_folders:
                avchains += chain +"|"
                avdsts = ""
                for dstver in checked_folders[chain]:
                    avdsts+=dstver+"|"
                    avconfigs = ""
                    for config in checked_folders[chain][dstver]:
                        avconfigs += config+"|"
                    pars.append(gammalib.GApplicationPar(chain+"_"+dstver+"_config","s","a",checked_folders[chain][dstver][0],avconfigs[:-1],"","Analysis configuration"))                
                pars.append(gammalib.GApplicationPar(chain+"_dstver","s","a",checked_folders[chain].keys()[0],avdsts[:-1],"","DST version"))       
            pars.append(gammalib.GApplicationPar("chain","s","a",checked_folders.keys()[0],avchains[:-1],"","Analysis chain"))
            pars.append(gammalib.GApplicationPar("runlistfile","f","a","runlist.lis","","","Runlist file"))   
            pars.append(gammalib.GApplicationPar("inmodel","f","h","NONE","","","Input model XML file (optional)"))
            pars.append(gammalib.GApplicationPar("outobs","f","a","obs.xml","","","Observation XML outfile"))
            pars.append(gammalib.GApplicationPar("outmodel","f","a","bgmodels.xml","","","Output model XML file"))
            pars.append(gammalib.GApplicationPar("datapath","s","h",self.datapath,"","","Data directoy"))     
            pars.append(gammalib.GApplicationPar("bkgpars","i","a","1","","","Number of free parameters per background model"))
            pars.append(gammalib.GApplicationPar("usetrig","b","h","yes","yes|no","","Pre-estimate background normalisation? Recommended for binned analysis"))
            pars.append(gammalib.GApplicationPar("bkgtype","s","a","irf","gauss|irf|aeff","","What background model should be used"))
            pars.append_standard()
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters       
        self.runlistfile = self["runlistfile"].filename()

        # Read models if user input
        self.inmodel = self["inmodel"].filename()
        if not (self.inmodel == "NONE" or self.inmodel == ""):
            self.inmodels = gammalib.GModels(self.inmodel)

        # Observation outfile
        self.outobs = self["outobs"].filename()

        # Output model file
        self.outmodel = self["outmodel"].filename()

        # Read hidden parameters
        self.bkgpars  = self["bkgpars"].integer()
        self.bkgtype  = self["bkgtype"].string()
        self.usetrig  = self["usetrig"].boolean()      
        self.datapath = self["datapath"].string()

        # Read run configuration parameters
        self.chain  = self["chain"].string()
        self.dstver = self[self.chain+"_dstver"].string()
        self.config = self[self.chain+"_"+self.dstver+"_config"].string()

        # Build directory name prefix
        self.directory = self.datapath+"/"
        self.directory += self.chain +"/"
        self.directory += self.dstver +"/"
        self.directory += self.config +"/"

        if not os.path.isdir(self.directory):
            raise gammalib.Exception.runtime_error("Folder structure not found. Please check your HESSFITS environment variable")

        self.m_log     = False # Logging in client tools
        self.m_debug   = False # Debugging in client tools
        self.m_clobber = self["clobber"].boolean()

        # Return
        return

    def check_structure(self,folders):

        # Initialise valid folders
        newfolders = {}

        # Check if we have found possible chains
        if len(folders):

            # Loop over chains
            for chain in folders:

                # Loop over dst versions
                for dst in folders[chain]:

                    # Loop over configs:
                    for config in folders[chain][dst]:

                        # check if config is valid
                        path = self.datapath+"/"+chain+"/"+dst+"/"+config
                        if os.path.isdir(path) and len(glob.glob(path+"/run*")):  
                            if not chain in newfolders:
                                newfolders[chain] = {}
                            if not dst in newfolders[chain]:
                                newfolders[chain][dst] = []
                            newfolders[chain][dst].append(config)

        # Return valid configs
        return newfolders
    
    def background_spectrum(self, run, prefactor, index, aefffile = ""):
        
        scale = pow(10,math.floor(math.log10(abs(prefactor))))        
        
        if index == 0.0 and self.bkgpars <= 1:
            spec = gammalib.GModelSpectralConst()
            spec["Value"].min(0.1*scale)
            spec["Value"].max(10*scale)
            spec["Value"].value(prefactor)
            spec["Value"].scale(scale)
            if self.bkgpars == 0:
                spec["Value"].fix()
            else:
                spec["Value"].free()
                
        else:
                
            # Create power law if number of free parameters are less 2 or less
            if self.bkgpars <= 2:
                e    = gammalib.GEnergy(1.0,"TeV")
                spec = gammalib.GModelSpectralPlaw(prefactor, index, e)  
                 
                # Set parameter ranges
                spec[0].scale(scale)
                spec[0].min(0.01 * scale)
                spec[0].max(100.0 * scale)
                spec[1].scale(1)
                spec[1].min(-5.0)
                spec[1].max(5.0)
                
                # Set number of free parameters
                if self.bkgpars == 0:
                    spec[0].fix()
                    spec[1].fix()
                elif self.bkgpars == 1:
                    spec[0].free()
                    spec[1].fix()
                else:
                    spec[0].free()
                    spec[1].free()           
            
            else:
                # Use several energy-dependent scaling factors  
                fits = gammalib.GFits(aefffile)
                emin = fits["EFFECTIVE AREA"].real("LO_THRES")
                emax = fits["EFFECTIVE AREA"].real("HI_THRES")
                fits.close()
                
                # Create reference powerlaw
                plaw = gammalib.GModelSpectralPlaw(prefactor, index, gammalib.GEnergy(1.0,"TeV")) 
                
                # Create spectral model and energy values
                spec = gammalib.GModelSpectralNodes()
                bounds = gammalib.GEbounds(self.bkgpars,gammalib.GEnergy(emin,"TeV"),gammalib.GEnergy(emax,"TeV"),True)
                for i in range(bounds.size()):     
                    energy = bounds.elogmean(i)
                    spec.append(energy, plaw.eval(energy, gammalib.GTime())) 
                for par in spec:
                    if "Energy" in par.name():
                        par.fix()
                    elif "Intensity" in par.name():  
                        par.scale(scale)                   
                        par.min(0.01*par.scale())
                        par.max(100.0*par.scale())
#         print prefactor, index
#         print spec
        return spec

    def hess_background(self, run, aefffile, trgrate, zenith):
        
        # Set IrfBackground
        if self.bkgtype == "irf":
            prefactor = 1.0
            index = 0.0
            if self.usetrig and trgrate > 0.0:
                # Use value from lookup
                prefactor = 0.00427233352611 * trgrate + 0.181244590076
            spec = self.background_spectrum(run, prefactor, index, aefffile)
            # Create background model instance
            bck = gammalib.GCTAModelIrfBackground(spec)
        
        # Set AeffBackground   
        elif self.bkgtype == "aeff":
            prefactor = 5e-14
            index = -2.5
            # Use values from lookup if possible
            if self.usetrig and trgrate > 0.0:
                prefactor = 1.52976410372e-16 * trgrate + 2.71805547108e-14               
            if zenith > 0.0: 
                index = 0.000271928721732 * zenith * zenith -0.0029602286453 * zenith -2.76404600416
        
            spec = self.background_spectrum(run, prefactor, index, aefffile)
                
            # Create background model instance
            bck = gammalib.GCTAModelAeffBackground(spec)
            
        # Set Gaussian Background
        elif self.bkgtype == "gauss":
            spec = self.background_spectrum(run, 1e-4, -1.8, aefffile)
            radial = gammalib.GCTAModelRadialGauss(2.5)
            bck = gammalib.GCTAModelRadialAcceptance(radial, spec)
        
        else:
            sys.exit("Background type \""+self.bkgtype+"\" unsupported")
        
        # Copy model
        model = bck.clone()
            
        # Assign specific run id
        model.ids(str(int(run)))
        
        # Assign instrument
        model.instruments("HESS")
        
        # Set name (arbitrary)
        model.name("bkg_"+run)  
        
        # Turn off TS calculation for background model
        model.tscalc(False)
        
        # Return model
        return model
        
    
    def erange(self):
        """
        Returns runlist energy range
        """
        return self.m_erange

    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()

        # Save residual map
        self.save(self.outobs, self.m_clobber)

        # Return
        return

    # Get runfolder from runnumber (subject to defined folder structure)
    def get_runfolder(self,runnr):

        start = runnr[:-2]+"00"
        if not int(runnr[-3]) % 2 == 0:
            start = str(int(start)-100)
            if len(start) == 5:
                start = "0"+start
        int_end = int(start)+199
        end = str(int_end)
        if len(end) == 5:
            end = "0"+ end

        runfold = runnr
        if len(runfold) == 5:
            runfold = "0"+runfold 

        return self.directory+"/run"+start+"-"+end+"/run"+runfold+"/"

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

        # Initialise empty models    
        self.models = gammalib.GModels()

        # Initialise empty xml file and append an observation list
        self.xml = gammalib.GXml()
        self.xml.append(gammalib.GXmlElement("observation_list title=\"observation list\""))
        lib = self.xml.element("observation_list", 0)

        # Read runlist from file
        runlist = []
        runfile = open(self.runlistfile)
        for line in runfile.readlines():
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if len(line.split()) > 0:
                runlist.append(line.split()[0])
            else:
                runlist.append('')
        runfile.close()

        # Loop over runs
        # Log output
        if self.logTerse():
            self.log("\n")
            self.log.header1("Looping over runs")

        # Initialise energy range values for logging
        runlist_emin = 100.0
        runlist_emax = 0.0

        for runnr in runlist:

            # Unify run number string length
            run = runnr
            if len(run) == 5:
                run = "0"+run

            # Get folder in which run files should be located
            runfolder = self.get_runfolder(run)

            # construct filenames
            eventfile = runfolder + "hess_events_"+run+".fits.gz"
            aefffile  = runfolder + "hess_aeff_2d_"+run+".fits.gz"
            psffile   = runfolder + "hess_psf_king_"+run+".fits.gz"
            edispfile = runfolder + "hess_edisp_2d_"+run+".fits.gz"
            bgfile    = runfolder + "hess_bkg_offruns_"+run+".fits.gz"

            # Check for existence of files
            msg = ""
            skip = False
            if not os.path.isdir(runfolder):
                msg = "Run "+str(int(run))+" does not exist in FITS format"
                skip = True
            elif not os.path.isfile(eventfile):
                msg = "Run "+str(int(run))+" has no eventfile - Run is skipped"
                skip = True
            elif not os.path.isfile(aefffile):
                msg = "Run "+str(int(run))+" has no effective area - Run is skipped"
                skip = True
            elif not os.path.isfile(psffile):
                msg = "Run "+str(int(run))+" has no PSF - Run is skipped"
                skip = True
            elif not os.path.isfile(edispfile):
                msg = "Run "+str(int(run))+" has no energy dispersion - Run is skipped"
                skip = True
            elif not os.path.isfile(bgfile):
                msg = "Run "+str(int(run))+" has no background - Run is skipped"
                skip = True            
            else:
                msg = "Adding run "+str(int(run))

            # Log ouput and skip run if not available         
            if self.logTerse():
                self.log(msg+"\n")
            if skip:
                continue
            if self.logExplicit():
                self.log(" Event file: "+eventfile+"\n")
                self.log(" Effective area: "+aefffile+"\n")
                self.log(" Point spread function: "+psffile+"\n")
                self.log(" Energy dispersion: "+edispfile+"\n")
                self.log(" Background file: "+bgfile+"\n")

            # Open fits file to determine the observation name
            fits = gammalib.GFits(eventfile)  
            object_name = fits["EVENTS"].string("OBJECT")

            # Retrieve trigger rate and zenith if available
            trgrate = 0.0
            zenith = 0.0
            if fits["EVENTS"].has_card("ZTRGRATE"):
                trgrate = fits["EVENTS"].real("ZTRGRATE")
            if fits["EVENTS"].has_card("ALT_PNT"):
                zenith = 90.-fits["EVENTS"].real("ALT_PNT")

            # Close FITS file
            fits.close()
            
            # Open effective area file to look for threshold
            aeff_fits = gammalib.GFits(aefffile)
            run_emin = aeff_fits["EFFECTIVE AREA"].real("LO_THRES")
            run_emax = aeff_fits["EFFECTIVE AREA"].real("HI_THRES")
            if run_emin < runlist_emin:
                runlist_emin = run_emin
            if run_emax > runlist_emax:
                runlist_emax = run_emax
            aeff_fits.close()          

            # Append observation to XML and set attributes
            obs = lib.append("observation");
            obs.attribute("name", object_name);
            obs.attribute("id", str(int(run)));
            obs.attribute("instrument", "HESS");

            # Append event file
            par = gammalib.GXmlElement("parameter name=\"EventList\"")
            par.attribute("file", eventfile);
            obs.append(par)

            # Append effective area
            aeff = gammalib.GXmlElement("parameter name=\"EffectiveArea\"")
            aeff.attribute("file",aefffile)

            # Append PSF
            psf = gammalib.GXmlElement("parameter name=\"PointSpreadFunction\"")
            psf.attribute("file",psffile)

            # Append energy dispersion
            edisp = gammalib.GXmlElement("parameter name=\"EnergyDispersion\"")
            edisp.attribute("file",edispfile)

            # Append background
            bck = gammalib.GXmlElement("parameter name=\"Background\"")
            bck.attribute("file",bgfile)

            # assign events and IRFs to observation
            obs.append(aeff)
            obs.append(psf)
            obs.append(edisp)
            obs.append(bck)
            
            # Append instrumental background model
            self.models.append(self.hess_background(run, aefffile, trgrate, zenith))

        # Continue only if there are observations available
        if lib.size():
            
            # Store energy range as member
            self.m_erange = gammalib.GEbounds(gammalib.GEnergy(runlist_emin,"TeV"),gammalib.GEnergy(runlist_emax,"TeV"))
               
            # Write energy range into log file    
            if self.logTerse():
                self.log("\n")
                self.log("Eenrgy range of obervation list: ")
                self.log(("% 3.2f" % runlist_emin)+" - "+("% 3.2f" % runlist_emax)+" TeV")
        
        else:
            self.log.header2("WARNING: No observation from given runlist available")

        # Append models provided by 'inmodels' if necessary
        if not self.inmodels == None:
            if self.logTerse():
                self.log("\n")
                self.log.header1("Appending models")

            # Loop over input models    
            for model in self.inmodels:
                if self.logTerse():
                    self.log.header3("Adding model \""+model.name()+"\"")
                self.models.append(model)

        # Return
        return

    def save(self,outfile,clobber):

        # Save observation XML file
        self.xml.save(outfile)

        # Save model XML file
        self.models.save(self.outmodel)


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates observation file.
    """
    # Create instance of application
    app = cshessobs(sys.argv)

    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
