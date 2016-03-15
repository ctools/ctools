#!/usr/bin/env python
# ==========================================================================
# Generation of an IACT observation definition file.
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

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name and version
        self._name    = "csiactobs"
        self._version = "1.1.0"

        # Initialise some members
        self._ebounds  = gammalib.GEbounds()
        self._datapath = os.getenv("VHEFITS","")
        self._inmodels = None
        self._xml      = gammalib.GXml()
        self._models   = gammalib.GModels()
        self._runlist  = []

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        
        if self._datapath == "":
            self._datapath = self["datapath"].string()
        
        # Check for input models
        if not self["inmodel"].filename() == "NONE":
            self._inmodels = gammalib.GModels(self["inmodel"].filename())
               
        # Expand environment
        self._datapath = gammalib.expand_env(self._datapath)
        
        # Read FITS production
        self._prodname    = self["prodname"].string()
        
        # Read runlist file if list not already filled
        if len(self._runlist) == 0:
            
            # Get file name
            self._runlistfile = self["infile"].filename()
            
            # Read runlist from file
            runfile = open(self._runlistfile.url())
            for line in runfile.readlines():
                if len(line) == 0:
                    continue
                if line[0] == '#':
                    continue
                if len(line.split()) > 0:
                    self._runlist.append(line.split()[0])
            runfile.close()
        
        # Read number of background parameters
        self._bkgpars     = self["bkgpars"].integer()
          
        # Output model file
        self._outmodel = self["outmodel"].filename()
        
        # Observation outfile
        self._outobs = self["outobs"].filename()
        
        # Master index file name
        self._master_indx = self["master_indx"].string()
        
        # Read flag for background scaling factor
        self._use_bkg_scale = self["bkg_scale"].boolean()
        
        # Read hierarchy of file loading
        self._ev_hiera      = self["ev_hiera"].string().split("|")
        self._aeff_hiera    = self["aeff_hiera"].string().split("|")
        self._psf_hiera     = self["psf_hiera"].string().split("|")
        self._bkg_hiera     = self["bkg_hiera"].string().split("|")  
        self._edisp_hiera   = self["edisp_hiera"].string().split("|") 
        self._bkg_mod_hiera = self["bkg_mod_hiera"].string().split("|")
        
        # Read hidden background parameters
        self._bkg_gauss_norm   = self["bkg_gauss_norm"].real()
        self._bkg_gauss_index  = self["bkg_gauss_index"].real()
        self._bkg_gauss_sigma  = self["bkg_gauss_sigma"].real()
        self._bkg_aeff_index   = self["bkg_aeff_index"].real()
        self._bkg_aeff_norm    = self["bkg_aeff_norm"].real()
        self._bkg_range_factor = self["bkg_range_factor"].real()
        
        # Open master index file and look for prodname 
        master_file = os.path.join(self._datapath, self._master_indx)
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
        self._hdu_index = self._obs_index = ""

        # Get HDUs
        for config in configs:
            if self._prodname == config["name"]:
                self._hdu_index = str(os.path.join(self._datapath, config["hduindx"]))
                self._obs_index = str(os.path.join(self._datapath, config["obsindx"]))
                break

        # Check HDUs
        if self._hdu_index == "" or self._obs_index == "":
            raise RuntimeError("*** ERROR: FITS data store \""+self._prodname+"\" not available. Run csiactdata to get a list of available storage names")
        filename = gammalib.GFilename(self._hdu_index+"[HDU_INDEX]")
        if not filename.is_fits():
            raise RuntimeError("*** ERROR: HDU index file \""+self._hdu_index+"[HDU_INDEX]\" for FITS data store \""+self._prodname+"\" not available. Check your master index file or run csiactdata to get a list of available storage names.")

        # Check for existence of "BKG_SCALE" in the observation index file if required
        if self._use_bkg_scale:
            filename = gammalib.GFilename(self._obs_index+"[OBS_INDEX]")
            if filename.is_fits():
                fits = gammalib.GFits(self._obs_index)
                if not fits["OBS_INDEX"].contains("BKG_SCALE"):
                    self._use_bkg_scale = False
                fits.close()
            else:
                self._use_bkg_scale = False
        
        # Create base data directory from hdu index file location
        self._subdir  = os.path.dirname(self._hdu_index)  
        self._debug   = False # Debugging in client tools
        self._clobber = self["clobber"].boolean()

        # Return
        return
    
    def runlist(self, runlist):
        """
        Set observation list
        """
        
        # Check if convertable to a python list
        try:
            
            # Store a copy of input argument
            self._runlist = list(runlist)
            
        except:
            
            # raise error if wrong type detected
            raise RuntimeError("Argument of csiactobs.runlist() must be of type 'list'")
        
        # Return
        return
    
    def _background_spectrum(self, run, prefactor, index, emin = 0.01, emax = 100.0):
        
        # Handle constant spectral model 
        if index == 0.0 and self._bkgpars <= 1:
            spec = gammalib.GModelSpectralConst()
            spec["Value"].min(prefactor / self._bkg_range_factor)
            spec["Value"].max(prefactor * self._bkg_range_factor)
            spec["Value"].value(prefactor)
            if self._bkgpars == 0:
                spec["Value"].fix()
            else:
                spec["Value"].free()
                
        else:
                
            # Create power law model 
            if self._bkgpars <= 2:
                e    = gammalib.GEnergy(1.0,"TeV")
                spec = gammalib.GModelSpectralPlaw(prefactor, index, e)  
                 
                # Set parameter ranges
                spec[0].min(0.01 * prefactor)
                spec[0].max(100.0 * prefactor)
                spec[1].scale(1)
                spec[1].min(-5.0)
                spec[1].max(5.0)
                
                # Set number of free parameters
                if self._bkgpars == 0:
                    spec[0].fix()
                    spec[1].fix()
                elif self._bkgpars == 1:
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
                bounds = gammalib.GEbounds(self._bkgpars,gammalib.GEnergy(emin,"TeV"),gammalib.GEnergy(emax,"TeV"), True)
                for i in range(bounds.size()):     
                    energy = bounds.elogmean(i)
                    value = plaw.eval(energy, gammalib.GTime())
                    spec.append(energy, value)
                for par in spec:
                    
                    if "Energy" in par.name():
                        par.fix()
                    elif "Intensity" in par.name():        
                        par.min(value / self._bkg_range_factor)
                        par.max(value * self._bkg_range_factor)
                        
        # Return spectrum
        return spec

    def _iact_background(self, telescope, obs_id, bkg_scale, bkgtype, emin=0.01, emax=100):
        
        # handle IrfBackground
        if bkgtype == "irf":
            prefactor  = 1.0
            index      = 0.0
            prefactor *= bkg_scale
            spec       = self._background_spectrum(obs_id, prefactor, index, emin, emax)
            
            # Create background model instance
            bck = gammalib.GCTAModelIrfBackground(spec)
        
        # Set AeffBackground   
        elif bkgtype == "aeff":
            prefactor = bkg_scale * self._bkg_aeff_norm
            spec      = self._background_spectrum(obs_id, prefactor, self._bkg_aeff_index, emin, emax)
                
            # Create background model instance
            bck = gammalib.GCTAModelAeffBackground(spec)
            
        # Set Gaussian Background
        elif bkgtype == "gauss":
            prefactor = bkg_scale * self._bkg_gauss_norm
            spec      = self._background_spectrum(obs_id, prefactor, self._bkg_gauss_index, emin, emax)
            radial    = gammalib.GCTAModelRadialGauss(self._bkg_gauss_sigma)
            bck       = gammalib.GCTAModelRadialAcceptance(radial, spec)
        
        else:
            sys.exit("Background type \""+self._bkgtype+"\" unsupported")
        
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

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Clear models
        self._models.clear()

        # Clear xml file and append an observation list
        self._xml.clear()
        self._xml.append(gammalib.GXmlElement("observation_list title=\"observation list\""))
        lib = self._xml.element("observation_list", 0)

        # Log output
        if self._logTerse():
            self._log("\n")
            self._log.header1("Looping over "+str(len(self._runlist))+" runs")

        # Initialise energy range values for logging
        self._ebounds.clear()

        # Loop over runs
        for obs_id in self._runlist:

            # Create selection string
            obs_selection = "[OBS_ID=="+str(obs_id)+"]"
            
            # Open HDU index file
            hduindx     = gammalib.GFits(self._hdu_index+"[HDU_INDEX]"+obs_selection)
            hduindx_hdu = hduindx["HDU_INDEX"]
            
            # Initialise files and hdu names
            eventfile = aefffile = psffile = edispfile = bkgfile = ""
            eventhdu  = aeffhdu  = ""
            
            types = []
            formats = []
            for i in range(hduindx_hdu.nrows()):  
                types.append(hduindx_hdu["HDU_TYPE"][i])
                formats.append(hduindx_hdu["HDU_CLASS"][i])
            
            # Handle events
            index = -1
            for version in self._ev_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                eventfile  = os.path.join(self._subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])
                eventhdu   = hduindx_hdu["HDU_NAME"][index]
                eventfile += "["+eventhdu+"]"
            if not gammalib.GFilename(eventfile).is_fits():
                if self._logTerse():
                    self._log("Skipping observation "+str(obs_id)+": eventfile \""+eventfile+"\" not found\n")
                continue
            
            # Handle Aeff
            index = -1
            for version in self._aeff_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                aefffile  = os.path.join(self._subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])
                aeffhdu   = hduindx_hdu["HDU_NAME"][index]
                aefffile += "["+aeffhdu+"]"
            if not gammalib.GFilename(aefffile).is_fits():
                if self._logTerse():
                    self._log("Skipping observation "+str(obs_id)+": effective area \""+aefffile+"\" not found\n")
                continue
                                
            # Handle psf
            index = -1
            for version in self._psf_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                psffile  = os.path.join(self._subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])
                psffile += "["+hduindx_hdu["HDU_NAME"][index]+"]"
            if not gammalib.GFilename(psffile).is_fits():
                if self._logTerse():
                    self._log("Skipping observation "+str(obs_id)+": point spread function \""+psffile+"\" not found\n")
                continue
                
            # Handle edisp
            index = -1
            for version in self._edisp_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                edispfile  = os.path.join(self._subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])   
                edispfile += "["+hduindx_hdu["HDU_NAME"][index]+"]"
            if not gammalib.GFilename(edispfile).is_fits():
                if self._logTerse():
                    self._log("Warning: observation "+str(obs_id)+" has no energy dispersion \""+edispfile+"\" information\n")
                    edispfile = ""
                
            # Handle background
            index = -1
            bkg_mod_hierarchy = list(self._bkg_mod_hiera)
            for version in self._bkg_hiera:
                n = formats.count(version) 
                if n > 0:
                    index = formats.index(version)
                    break
            if index >= 0:
                bkgfile  = os.path.join(self._subdir, hduindx_hdu["FILE_DIR"][index], hduindx_hdu["FILE_NAME"][index])    
                bkgfile += "["+hduindx_hdu["HDU_NAME"][index]+"]"
            if not gammalib.GFilename(bkgfile).is_fits():
                if self._logTerse():
                    self._log("Warning: observation "+str(obs_id)+" has no background information (file=\""+bkgfile+"\"). IRF background cannot be used\n")
                    bkgfile = ""
                    if "irf" in bkg_mod_hierarchy:
                        bkg_mod_hierarchy.remove("irf")

            # Close hdu index file
            hduindx.close()

            # Handle background scale information if available
            bkg_scale = 1.0
            if self._use_bkg_scale:
                obsindx = gammalib.GFits(self._obs_index+"[OBS_INDEX]"+obs_selection)
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
            aeff_fits  = gammalib.GFits(aefffile)
            aeff_table = aeff_fits[aeffhdu]
            
            # Set energy range from header keyword if present
            if aeff_table.has_card("LO_THRES") and aeff_table.has_card("HI_THRES"):
                
                # Get safe energy range
                run_emin = gammalib.GEnergy(aeff_table.real("LO_THRES"), "TeV")
                run_emax = gammalib.GEnergy(aeff_table.real("HI_THRES"), "TeV")
                
                # Append to ebounds
                self._ebounds.append(run_emin, run_emax)
            
            else:
                # Set default values for energy range
                run_emin = gammalib.GEnergy(10,"GeV")
                run_emax = gammalib.GEnergy(100,"TeV")
            
            # Close Aeff fits file
            aeff_fits.close()     
            
            # Append instrumental background model
            self._models.append(self._iact_background(telescope, obs_id, bkg_scale, bkg_mod_hierarchy[0], run_emin.TeV(), run_emax.TeV()))

            # Logging
            if self._logTerse():
                self._log("Adding observation "+str(obs_id)+" (\""+object_name+"\")\n")

            if self._logExplicit():
                self._log(" Event file: "+eventfile+"\n")
                self._log(" Effective area: "+aefffile+"\n")
                self._log(" Point spread function: "+psffile+"\n")
                self._log(" Energy dispersion: "+edispfile+"\n")
                self._log(" Background file: "+bkgfile+"\n")  
                if self._use_bkg_scale:
                    self._log(" Background scale: "+str(bkg_scale)+"\n")  
                self._log("\n")

            # Append observation to XML and set attributes
            obs = lib.append("observation")
            obs.attribute("name", object_name)
            obs.attribute("id", str(obs_id))
            obs.attribute("instrument", telescope)

            # Append event file
            ev = gammalib.GXmlElement("parameter name=\"EventList\"")
            ev.attribute("file", eventfile)

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
            if self._logTerse():
                self._log("\n")
                self._log.header3("Energy range of obervation list")
            
            # Logging if energy range is available
            if self._ebounds.size() > 0:
                               
                # Write energy range into log file    
                self._log(str(self._ebounds.emin())+" - "+str(self._ebounds.emax()))
                self._log("\n")
                
            else:
               
                # Write 'not available' into log file    
                self._log("not available")
                self._log("\n")
        
        else:
            self._log.header2("WARNING: No observation from given runlist available")

        # Append models provided by 'inmodels' if necessary
        if not self._inmodels == None:
            if self._logTerse():
                self._log("\n")
                self._log.header1("Appending models")

            # Loop over input models    
            for model in self._inmodels:
                if self._logTerse():
                    self._log.header3("Adding model \""+model.name()+"\"")
                self._models.append(model)

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

        # Save residual map
        self.save(self._outobs, self._clobber)

        # Return
        return

    def save(self,outfile,clobber):

        # Save observation XML file
        self._xml.save(outfile)

        # Save model XML file
        self._models.save(self._outmodel)

        # Return
        return       

    def ebounds(self):
        """
        Returns runlist energy range
        """
        return self._ebounds
    
    def obs(self):
        """
        Returns GObservations object
        """   
        # Initialise observations
        obs = gammalib.GObservations()
        
        # Read observations if XML is filled
        if self._xml.size():   
            obs.read(self._xml)
                
        # Assign models
        obs.models(self._models)
        
        # Return observations
        return obs


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csiactobs(sys.argv)

    # Execute application
    app.execute()
