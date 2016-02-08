#!/usr/bin/env python
# ==========================================================================
# Generation of an IACT observation definition file.
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


# =============== #
# csiactobs class #
# =============== #
class csiactobs(ctools.cscript):
    """
    This class implements the creation of a observation xml file for IACT
    data analysis. This class is dedicated for use inside a IACT
    Collaboration, i.e. it can only be used if you have access to IACT data
    in FITS format. The FITS data has to be structured and in the format
    described here:
    http://gamma-astro-data-formats.readthedocs.org/en/latest/
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name and version
        self.name    = "csiactobs"
        self.version = "1.1.0"

        # Initialise some members
        self.m_ebounds = gammalib.GEbounds()
        self.datapath = os.getenv("VHEFITS","")
        self.inmodels = None
        self.outobs   = "obs.xml"

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
            # Signal if parfile was not found
            sys.stdout.write("Parfile "+parfile+" not found. Create default parfile.\n")

            # Create default parfile
            pars = gammalib.GApplicationPars()    
            pars.append(gammalib.GApplicationPar("datapath","s","a",self.datapath,"","","Path were data is located"))       
            pars.append(gammalib.GApplicationPar("prodname","s","a","","","","Data storage name (Run csiactdata to view your options)"))
            pars.append(gammalib.GApplicationPar("infile","f","a","runlist.lis","","","Runlist file"))   
            pars.append(gammalib.GApplicationPar("inmodel","f","h","NONE","","","Input model XML file (optional)"))
            pars.append(gammalib.GApplicationPar("outobs","f","a","obs.xml","","","Observation XML outfile"))
            pars.append(gammalib.GApplicationPar("outmodel","f","a","bgmodels.xml","","","Output model XML file"))
            pars.append(gammalib.GApplicationPar("bkgpars","i","a","1","","","Number of free parameters per background model"))
            pars.append(gammalib.GApplicationPar("master_indx","s","h","master.json","","","Name of master index file"))
            pars.append(gammalib.GApplicationPar("bkg_scale","b","h","yes","","","Specifies whether the background scaling factor from the observation index file should be applied if available. "))
            pars.append(gammalib.GApplicationPar("ev_hiera","s","h","events","","","Hierarchy of event formats"))
            pars.append(gammalib.GApplicationPar("aeff_hiera","s","h","aeff_2d","","","Hierarchy of effective area formats"))
            pars.append(gammalib.GApplicationPar("psf_hiera","s","h","psf_king|psf_3gauss","","","Hierarchy of psf formats"))
            pars.append(gammalib.GApplicationPar("edisp_hiera","s","h","edisp_2d","","","Hierarchy of energy dispersion formats"))
            pars.append(gammalib.GApplicationPar("bkg_hiera","s","h","bkg_3d","","","Hierarchy of background formats"))
            pars.append(gammalib.GApplicationPar("bkg_mod_hiera","s","h","irf|aeff|gauss","","","Hierarchy of background models"))
            pars.append(gammalib.GApplicationPar("bkg_gauss_norm","r","h","1e-8","","","Input normalisation for Gaussian background"))
            pars.append(gammalib.GApplicationPar("bkg_gauss_index","r","h","-2.0","","","Input spectral index for Gaussian background"))
            pars.append(gammalib.GApplicationPar("bkg_gauss_sigma","r","h","2.5","","","Input sigma for Gaussian background"))
            pars.append(gammalib.GApplicationPar("bkg_aeff_norm","r","h","1e-14","","","Input normalisation for effective area background"))
            pars.append(gammalib.GApplicationPar("bkg_aeff_index","r","h","-2.0","","","Input spectral index for effective area background"))
            pars.append(gammalib.GApplicationPar("bkg_range_factor","r","h","100.0","","","Factor to determine range of background normalisation"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csiactobs.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        
        if self.datapath == "":
            self.datapath = self["datapath"].string()
        
        # Check for input models
        if not self["inmodel"].filename() == "NONE":
            self.inmodels = gammalib.GModels(self["inmodel"].filename())
               
        # Expand environment
        self.datapath = gammalib.expand_env(self.datapath)
        
        # Read user parameters  
        self.m_prodname    = self["prodname"].string()
        self.m_runlistfile = gammalib.expand_env(self["infile"].filename())
        self.m_bkgpars     = self["bkgpars"].integer()
          
        # Output model file
        self.outmodel = self["outmodel"].filename()
        
        # Observation outfile
        self.outobs = self["outobs"].filename()
        
        # Master index file name
        self.m_master_indx = self["master_indx"].string()
        
        # Read flag for background scaling factor
        self.m_use_bkg_scale = self["bkg_scale"].boolean()
        
        # Read hierarchy of file loading
        self.m_ev_hiera      = self["ev_hiera"].string().split("|")
        self.m_aeff_hiera    = self["aeff_hiera"].string().split("|")
        self.m_psf_hiera     = self["psf_hiera"].string().split("|")
        self.m_bkg_hiera     = self["bkg_hiera"].string().split("|")  
        self.m_edisp_hiera   = self["edisp_hiera"].string().split("|") 
        self.m_bkg_mod_hiera = self["bkg_mod_hiera"].string().split("|")
        
        # Read hidden background parameters
        self.m_bkg_gauss_norm   = self["bkg_gauss_norm"].real()
        self.m_bkg_gauss_index  = self["bkg_gauss_index"].real()
        self.m_bkg_gauss_sigma  = self["bkg_gauss_sigma"].real()
        self.m_bkg_aeff_index   = self["bkg_aeff_index"].real()
        self.m_bkg_aeff_norm    = self["bkg_aeff_norm"].real()
        self.m_bkg_range_factor = self["bkg_range_factor"].real()
        
        # Open master index file and look for prodname 
        master_file = os.path.join(self.datapath, self.m_master_indx)
        if not os.path.isfile(master_file):
            raise RuntimeError("FITS data store not available. No master index file found. Make sure the file is copied from the server and your datapath is set correctly.")

        # Open and load JSON file
        json_data = open(master_file).read()
        data      = json.loads(json_data)    
        if not "datasets" in data:
            raise RuntimeError("Key \"datasets\" not available in master index file.")

        # Get configurations
        configs = data["datasets"]

        # Initialise HDUs
        self.m_hdu_index = self.m_obs_index = ""

        # Get HDUs
        for config in configs:
            if self.m_prodname == config["name"]:
                self.m_hdu_index = str(os.path.join(self.datapath, config["hduindx"]))
                self.m_obs_index = str(os.path.join(self.datapath, config["obsindx"]))
                break

        # Check HDUs
        if self.m_hdu_index == "" or self.m_obs_index == "":
            raise RuntimeError("*** ERROR: FITS data store \""+self.m_prodname+"\" not available. Run csiactdata to get a list of available storage names")
        if not gammalib.is_fits(self.m_hdu_index+"[HDU_INDEX]"):
            raise RuntimeError("*** ERROR: HDU index file \""+self.m_hdu_index+"[HDU_INDEX]\" for FITS data store \""+self.m_prodname+"\" not available. Check your master index file or run csiactdata to get a list of available storage names.")

        # Check for existence of "BKG_SCALE" in the observation index file if required
        if self.m_use_bkg_scale:
            if gammalib.is_fits(self.m_obs_index+"[OBS_INDEX]"):
                fits = gammalib.GFits(self.m_obs_index)
                if not fits["OBS_INDEX"].contains("BKG_SCALE"):
                    self.m_use_bkg_scale = False
                fits.close()
            else:
                self.m_use_bkg_scale = False
        
        # Create base data directory from hdu index file location
        self.subdir    = os.path.dirname(self.m_hdu_index)  
        self.m_log     = False # Logging in client tools
        self.m_debug   = False # Debugging in client tools
        self.m_clobber = self["clobber"].boolean()

        # Return
        return
    
    def background_spectrum(self, run, prefactor, index, emin = 0.01, emax = 100.0):
        
        # Handle constant spectral model 
        if index == 0.0 and self.m_bkgpars <= 1:
            spec = gammalib.GModelSpectralConst()
            spec["Value"].min(prefactor / self.m_bkg_range_factor)
            spec["Value"].max(prefactor * self.m_bkg_range_factor)
            spec["Value"].value(prefactor)
            if self.m_bkgpars == 0:
                spec["Value"].fix()
            else:
                spec["Value"].free()
                
        else:
                
            # Create power law model 
            if self.m_bkgpars <= 2:
                e    = gammalib.GEnergy(1.0,"TeV")
                spec = gammalib.GModelSpectralPlaw(prefactor, index, e)  
                 
                # Set parameter ranges
                spec[0].min(0.01 * prefactor)
                spec[0].max(100.0 * prefactor)
                spec[1].scale(1)
                spec[1].min(-5.0)
                spec[1].max(5.0)
                
                # Set number of free parameters
                if self.m_bkgpars == 0:
                    spec[0].fix()
                    spec[1].fix()
                elif self.m_bkgpars == 1:
                    spec[0].free()
                    spec[1].fix()
                else:
                    spec[0].free()
                    spec[1].free()           
            
            else:
                
                # Create reference powerlaw
                plaw = gammalib.GModelSpectralPlaw(prefactor, index, gammalib.GEnergy(1.0,"TeV")) 
                
                # Create spectral model and energy values
                spec = gammalib.GModelSpectralNodes()
                bounds = gammalib.GEbounds(self.m_bkgpars,gammalib.GEnergy(emin,"TeV"),gammalib.GEnergy(emax,"TeV"), True)
                for i in range(bounds.size()):     
                    energy = bounds.elogmean(i)
                    value = plaw.eval(energy, gammalib.GTime())
                    spec.append(energy, value)
                for par in spec:
                    
                    if "Energy" in par.name():
                        par.fix()
                    elif "Intensity" in par.name():        
                        par.min(value / self.m_bkg_range_factor)
                        par.max(value * self.m_bkg_range_factor)
                        
        # Return spectrum
        return spec

    def iact_background(self, telescope, obs_id, bkg_scale, bkgtype, emin=0.01, emax=100):
        
        # handle IrfBackground
        if bkgtype == "irf":
            prefactor  = 1.0
            index      = 0.0
            prefactor *= bkg_scale
            spec = self.background_spectrum(obs_id, prefactor, index, emin, emax)
            
            # Create background model instance
            bck = gammalib.GCTAModelIrfBackground(spec)
        
        # Set AeffBackground   
        elif bkgtype == "aeff":
        
            prefactor = bkg_scale * self.m_bkg_aeff_norm
            spec = self.background_spectrum(obs_id, prefactor, self.m_bkg_aeff_index, emin, emax)
                
            # Create background model instance
            bck = gammalib.GCTAModelAeffBackground(spec)
            
        # Set Gaussian Background
        elif bkgtype == "gauss":
            prefactor = bkg_scale * self.m_bkg_gauss_norm
            spec      = self.background_spectrum(obs_id, prefactor, self.m_bkg_gauss_index, emin, emax)
            radial    = gammalib.GCTAModelRadialGauss(self.m_bkg_gauss_sigma)
            bck       = gammalib.GCTAModelRadialAcceptance(radial, spec)
        
        else:
            sys.exit("Background type \""+self.bkgtype+"\" unsupported")
        
        # Copy model
        model = bck.clone()
            
        # Assign specific run id
        model.ids(str(obs_id))
        
        # Assign instrument
        model.instruments(telescope)
        
        # Set name (arbitrary)
        model.name("bkg_"+str(obs_id))  
        
        # Turn off TS calculation for background model
        model.tscalc(False)
        
        # Return model
        return model

    def ebounds(self):
        """
        Returns runlist energy range
        """
        return self.m_ebounds

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
        runfile = open(self.m_runlistfile)
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

        # Log output
        if self.logTerse():
            self.log("\n")
            self.log.header1("Looping over "+str(len(runlist))+" runs")

        # Initialise energy range values for logging
        self.m_ebounds.clear()

        # Loop over runs
        for obs_id in runlist:

            # Create selection string
            obs_selection = "[OBS_ID=="+str(obs_id)+"]"
            
            # Open HDU index file
            hduindx = gammalib.GFits(self.m_hdu_index+"[HDU_INDEX]"+obs_selection)
            hduindx_hdu = hduindx["HDU_INDEX"]
            
            # Initialise files and hdu names
            eventfile = aefffile = psffile = edispfile = bkgfile = ""
            eventhdu = aeffhdu = ""
            
            types = []
            formats = []
            for i in range(hduindx_hdu.nrows()):  
                types.append(hduindx_hdu["HDU_TYPE"][i])
                formats.append(hduindx_hdu["HDU_CLASS"][i])
            
            # Handle events
            index = -1
            for version in self.m_ev_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                eventfile = os.path.join(self.subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])
                eventhdu = hduindx_hdu["HDU_NAME"][index]
                eventfile += "["+eventhdu+"]"
            if not gammalib.is_fits(eventfile):
                if self.logTerse():
                    self.log("Skipping observation "+str(obs_id)+": eventfile \""+eventfile+"\" not found\n")
                continue
            
            # Handle Aeff
            index = -1
            for version in self.m_aeff_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                aefffile = os.path.join(self.subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])
                aeffhdu  = hduindx_hdu["HDU_NAME"][index]
                aefffile += "["+aeffhdu+"]"
            if not gammalib.is_fits(aefffile):
                if self.logTerse():
                    self.log("Skipping observation "+str(obs_id)+": effective area \""+aefffile+"\" not found\n")
                continue
                                
            # Handle psf
            index = -1
            for version in self.m_psf_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                psffile = os.path.join(self.subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])
                psffile += "["+hduindx_hdu["HDU_NAME"][index]+"]"
            if not gammalib.is_fits(psffile):
                if self.logTerse():
                    self.log("Skipping observation "+str(obs_id)+": point spread function \""+psffile+"\" not found\n")
                continue
                
            # Handle edisp
            index = -1
            for version in self.m_edisp_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                edispfile = os.path.join(self.subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])   
                edispfile += "["+hduindx_hdu["HDU_NAME"][index]+"]"
            if not gammalib.is_fits(edispfile):
                if self.logTerse():
                    self.log("Warning: observation "+str(obs_id)+" has no energy dispersion \""+edispfile+"\" information\n")
                    edispfile = ""
                
            # Handle background
            index = -1
            bkg_mod_hierarchy = self.m_bkg_mod_hiera
            for version in self.m_bkg_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                bkgfile = os.path.join(self.subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])    
                bkgfile += "["+hduindx_hdu["HDU_NAME"][index]+"]"
            if not gammalib.is_fits(bkgfile):
                if self.logTerse():
                    self.log("Warning: observation "+str(obs_id)+" has no background information (file=\""+bkgfile+"\"). IRF background cannot be used\n")
                    bkgfile = ""
                    bkg_mod_hierarchy.remove("irf")

            # Close hdu index file
            hduindx.close()

            # Handle background scale information if available
            bkg_scale = 1.0
            if self.m_use_bkg_scale:
                obsindx = gammalib.GFits(self.m_obs_index+"[OBS_INDEX]"+obs_selection)
                bkg_scale = obsindx["OBS_INDEX"]["BKG_SCALE"][0]
                obsindx.close()
            
            # Open fits file to determine the observation name
            fits = gammalib.GFits(eventfile)  
            events = fits[eventhdu]      
            object_name = events.string("OBJECT")
            telescope = events.string("TELESCOP")

            # Close FITS file
            fits.close()
            
            # Open effective area file to look for threshold
            aeff_fits = gammalib.GFits(aefffile)
            aeff_table = aeff_fits[aeffhdu]
            
            # Set energy range from header keyword if present
            if aeff_table.has_card("LO_THRES") and aeff_table.has_card("HI_THRES"):
                
                # Get safe energy range
                run_emin = gammalib.GEnergy(aeff_table.real("LO_THRES"), "TeV")
                run_emax = gammalib.GEnergy(aeff_table.real("HI_THRES"), "TeV")
                
                # Append to ebounds
                self.m_ebounds.append(run_emin, run_emax)
            
            else:
                # Set default values for energy range
                run_emin = gammalib.GEnergy(10,"GeV")
                run_emax = gammalib.GEnergy(100,"TeV")
            
            # Close Aeff fits file
            aeff_fits.close()     
            
            # Append instrumental background model
            self.models.append(self.iact_background(telescope, obs_id, bkg_scale, bkg_mod_hierarchy[0], run_emin.TeV(), run_emax.TeV()))

            # Logging
            if self.logTerse():
                self.log("Adding observation "+str(obs_id)+" (\""+object_name+"\")\n")

            if self.logExplicit():
                self.log(" Event file: "+eventfile+"\n")
                self.log(" Effective area: "+aefffile+"\n")
                self.log(" Point spread function: "+psffile+"\n")
                self.log(" Energy dispersion: "+edispfile+"\n")
                self.log(" Background file: "+bkgfile+"\n")  
                if self.m_use_bkg_scale:
                    self.log(" Background scale: "+str(bkg_scale)+"\n")  
                self.log("\n")

            # Append observation to XML and set attributes
            obs = lib.append("observation");
            obs.attribute("name", object_name);
            obs.attribute("id", str(obs_id));
            obs.attribute("instrument", telescope);

            # Append event file
            ev = gammalib.GXmlElement("parameter name=\"EventList\"")
            ev.attribute("file", eventfile);

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
            bck.attribute("file",bkgfile)

            # assign events and IRFs to observation
            obs.append(ev)
            obs.append(aeff)
            obs.append(psf)
            obs.append(edisp)
            obs.append(bck)
            
        # Continue only if there are observations available
        if lib.size():
            
            # Log header of energy range
            if self.logTerse():
                self.log("\n")
                self.log.header3("Eenrgy range of obervation list: ")
            
            # Logging if energy range is available
            if self.m_ebounds.size() > 0:
                               
                # Write energy range into log file    
                self.log(str(self.m_ebounds.emin())+" - "+str(self.m_ebounds.emax()))
                self.log("\n")
                
            else:
               
                # Write 'not available' into log file    
                self.log("not available")
                self.log("\n")
        
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

        # Return
        return       


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates observation file.
    """
    # Create instance of application
    app = csiactobs(sys.argv)

    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
