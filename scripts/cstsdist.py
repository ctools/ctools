#! /usr/bin/env python
# ==========================================================================
# This script generates the TS distribution for a particular model based
# on Monte-Carlo simulations.
#
# Copyright (C) 2011-2015 Juergen Knoedlseder
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
from ctools import obsutils
import sys
import csv


# ============== #
# cstsdist class #
# ============== #
class cstsdist(ctools.cscript):
    """
    This class implements the TS distribution generation script. It derives
    from the GammaLib::GApplication class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "cstsdist"
        self.version = "0.2.0"
        
        # Initialise some members
        self.obs        = None
        self.bkg_model  = None
        self.full_model = None
        
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
            pars.append(gammalib.GApplicationPar("inobs","f","h","NONE","","","Event list, counts cube, or observation definition file"))
            pars.append(gammalib.GApplicationPar("outfile","f","h","ts.dat","","","Output file name"))
            pars.append(gammalib.GApplicationPar("expcube","s","a","NONE","","","Exposure cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("psfcube","s","a","NONE","","","PSF cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","dummy","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","cta_dummy_irf","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("ra","r","h","83.6331","0","360","RA of pointing (deg)"))
            pars.append(gammalib.GApplicationPar("dec","r","h","22.0145","-90","90","Dec of pointing (deg)"))
            pars.append(gammalib.GApplicationPar("emin","r","a","0.1","","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","a","100.0","","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("tmin","r","h","0.0","","","Start time (MET in s)"))
            pars.append(gammalib.GApplicationPar("tmax","r","a","1800.0","","","Observation duration (in s)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","0","","","Number of energy bins (0=unbinned)"))
            pars.append(gammalib.GApplicationPar("npix","i","h","200","","","Number of pixels for binned"))
            pars.append(gammalib.GApplicationPar("binsz","r","h","0.05","","","Pixel size for binned (deg/pixel)"))
            pars.append(gammalib.GApplicationPar("deadc","r","h","0.95","","","Deadtime correction factor"))
            pars.append(gammalib.GApplicationPar("rad","r","h","5.0","","","Radius of ROI (deg)"))         
            pars.append(gammalib.GApplicationPar("ntrials","i","a","10","","","Number of trials"))
            pars.append(gammalib.GApplicationPar("type","s","a","point","","","Source model type (point/gauss/shell/disk)"))
            pars.append(gammalib.GApplicationPar("index","r","h","-2.48","","","Spectral index"))       
            pars.append(gammalib.GApplicationPar("offset","r","a","0.0","0.0","","Source offset angle (deg)"))
            pars.append(gammalib.GApplicationPar("bkg","s","a","$CTOOLS/share/models/bkg_dummy.txt","","","Background model file function"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","cstsdist.log","","","Log filename"))
            pars.save(parfile)
        
        # Return
        return
        
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        
        # Set observation if not done before
        if self.obs == None or self.obs.size() == 0:
            self.obs = self.get_observations()
            
        # Get number of energy bins
        self.m_enumbins = self["enumbins"].integer()
        
        # Read parameters for binned if requested
        if not self.m_enumbins == 0:
            self.m_npix     = self["npix"].integer()
            self.m_binsz    = self["binsz"].real()
        else:
            # Set dummy values (required by obsutils)
            self.m_npix = 0
            self.m_binsz = 0.0
          
        # Get other parameters
        self.m_outfile  = self["outfile"].filename()
        self.m_ntrials  = self["ntrials"].integer()
        self.m_type     = self["type"].string()
        self.m_index    = self["index"].real()
        self.m_offset   = self["offset"].real()
        self.m_bkg      = self["bkg"].string()

        # Set some fixed parameters
        self.m_log   = False # Logging in client tools
        self.m_debug = False # Debugging in client tools
          
        pnt = gammalib.GSkyDir()
        pnt.radec_deg(self["ra"].real(),self["dec"].real())  
          
        # Initialise models. Note that we centre the point source at the
        # center of our observation, so we're onaxis.
        self.bkg_model  = gammalib.GModels()
        self.full_model = gammalib.GModels()
        self.bkg_model.append(self.set_bkg_model())
        self.full_model.append(self.set_bkg_model())
        self.full_model.append(self.set_src_model(pnt.l_deg(), pnt.b_deg()+self.m_offset, \
                                                  flux=0.010, \
                                                  type=self.m_type, \
                                                  index=self.m_index))

        # Attach background model to observation container
        self.obs.models(self.bkg_model)
        
        # Return
        return
    
    def models(self, models):
        """
        Set model.
        """
        # Copy models
        self.obs.models(models.copy())
    
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

        # Write models into logger
        if self.logTerse():
            self.log("\n")
            self.log.header1("Test model")
            self.log(str(self.full_model))
            self.log("\n")

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Generate TS distribution")

        # Loop over trials
        for seed in range(self.m_ntrials):
        
            # Make a trial
            result = self.trial(seed)
            
            # Write out result immediately
            if seed == 0:
                file = open(self.m_outfile, 'w')
                writer = csv.DictWriter(file, result['colnames'])
                writer.writerow(dict((_,_) for _ in result['colnames']))
            else:
                file = open(self.m_outfile, 'a')
            writer = csv.DictWriter(file, result['colnames'])
            writer.writerow(result['values'])
            file.close()
        
        # Return
        return
  
    def set_bkg_model(self, fitsigma=False):
        """
        Setup CTA background model.
        """
        # Define radial component
        radial = gammalib.GCTAModelRadialGauss(3.0)
        if fitsigma:
            radial["Sigma"].free()
        else:
            radial["Sigma"].fix()
        
        # Define spectral component
        spectrum = gammalib.GModelSpectralFunc(self.m_bkg, 1.0)
        
        # Create background model
        model = gammalib.GCTAModelRadialAcceptance(radial, spectrum)
        model.name("Background")
        model.instruments("CTA")
    
        # Return background model
        return model
    
    def set_src_model(self, l, b, flux=1.0, index=-2.48, \
                      type="point", sigma=1.0, radius=1.0, width=0.1, \
                      fitpos=False, fitidx=False):
        """
        Returns a single source with Crab-like spectrum. The source flux
        can be scaled in Crab units. The Crab spectrum is based on MAGIC
        observations (Albert et al. 2008, ApJ, 674, 1037).

        Parameters:
         l      - Galactic longitude of source location [deg]
         b      - Galactic latitude of source location [deg]
        Keywords:
         flux   - Source flux [Crabs]
         index  - Spectral index
         type   - Source type ("point", "gauss", "disk", "shell")
         sigma  - Gaussian sigma (for type="gauss")
         radius - Disk or shell inner radius [deg] (for type="disk" and type="shell")
         width  - Shell width [deg] (for type="shell")
         fitpos - Fit position? (default: True)
         fitidx - Fit index? (default: True)
        """
        # Set source location
        location = gammalib.GSkyDir()
        location.lb_deg(l, b)
    
        # Set source spectrum
        spectrum = gammalib.GModelSpectralPlaw(flux*5.7e-16, index, gammalib.GEnergy(0.3, "TeV"))
        if fitidx:
            spectrum["Index"].free()
        else:
            spectrum["Index"].fix() 

        # Set source
        if type == "point":
            spatial = gammalib.GModelSpatialPointSource(location)
            if fitpos:
                spatial[0].free()
                spatial[1].free()
        elif type == "gauss":
            spatial = gammalib.GModelSpatialRadialGauss(location, sigma)
            if fitpos:
                spatial[0].free()
                spatial[1].free()
        elif type == "disk":
            spatial = gammalib.GModelSpatialRadialDisk(location, radius)
            if fitpos:
                spatial[0].free()
                spatial[1].free()
        elif type == "shell":
            spatial = gammalib.GModelSpatialRadialShell(location, radius, width)
            if fitpos:
                spatial[0].free()
                spatial[1].free()
        else:
            self.log("ERROR: Unknown source type '"+type+"'.\n")
            return None
        source = gammalib.GModelSky(spatial, spectrum)

        # Set source name
        source.name("Test")
        
        # Turn TS value computation off
        # since we compute it by hand later on
        source.tscalc(False)
    
        # Return source
        return source
    
    def trial(self, seed):
        """
        Create the TS for a single trial.
        
        Parameters:
         seed - Random number generator seed
        """
        # Write header
        if self.logExplicit():
            self.log.header2("Trial "+str(seed+1))

        # Simulate events
        sim = obsutils.sim(self.obs, \
                           nbins=self.m_enumbins, \
                           seed=seed, \
                           binsz=self.m_binsz, \
                           npix=self.m_npix, \
                           log=self.m_log, debug=self.m_debug)
        
        # Determine number of events in simulation
        nevents = 0.0
        for run in sim:
            nevents += run.events().number()

        # Write simulation results
        if self.logExplicit():
            self.log.header3("Simulation")
            self.log.parformat("Number of simulated events")
            self.log(nevents)
            self.log("\n")
        
        # Fit background only
        sim.models(self.bkg_model)
        like_bgm   = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
        result_bgm = like_bgm.obs().models()
        LogL_bgm   = like_bgm.opt().value()
        npred_bgm  = like_bgm.obs().npred()

        # Write background fit results
        if self.logExplicit():
            self.log.header3("Background model fit")
            self.log.parformat("log likelihood")
            self.log(LogL_bgm)
            self.log("\n")
            self.log.parformat("Number of predicted events")
            self.log(npred_bgm)
            self.log("\n")
            for model in result_bgm:
                self.log.parformat("Model")
                self.log(model.name())
                self.log("\n")
                for par in model:
                    self.log(str(par)+"\n")

        # Fit background and test source
        sim.models(self.full_model)
        like_all   = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
        result_all = like_all.obs().models()
        LogL_all   = like_all.opt().value()
        npred_all  = like_all.obs().npred()
        ts         = 2.0*(LogL_bgm-LogL_all)
        
        # Write background and test source fit results
        if self.logExplicit():
            self.log.header3("Background and test source model fit")
            self.log.parformat("Test statistics")
            self.log(ts)
            self.log("\n")
            self.log.parformat("log likelihood")
            self.log(LogL_all)
            self.log("\n")
            self.log.parformat("Number of predicted events")
            self.log(npred_all)
            self.log("\n")
            for model in result_all:
                self.log.parformat("Model")
                self.log(model.name())
                self.log("\n")
                for par in model:
                    self.log(str(par)+"\n")
                    
        # Write result
        elif self.logTerse():
            self.log.parformat("Trial "+str(seed))
            self.log("TS=")
            self.log(ts)
            self.log("  Prefactor=")
            self.log(result_all["Test"]["Prefactor"].value())
            self.log("+/-")
            self.log(result_all["Test"]["Prefactor"].error())
            self.log("\n")
        
        # Initialise results
        colnames = []
        values   = {}
        
        # Set TS value
        colnames.append("TS")
        values["TS"] = ts

        # Set logL for background fit
        colnames.append("LogL_bgm")
        values["LogL_bgm"] = LogL_bgm

        # Set logL for full fit
        colnames.append("LogL_all")
        values["LogL_all"] = LogL_all

        # Set Nevents
        colnames.append("Nevents")
        values["Nevents"] = nevents

        # Set Npred for background fit
        colnames.append("Npred_bkg")
        values["Npred_bkg"] = npred_bgm

        # Set Npred for full fit
        colnames.append("Npred_all")
        values["Npred_all"] = npred_all
        
        # Gather free full fit parameters
        for i in range(result_all.size()):
            model      = result_all[i]
            model_name = model.name()
            for k in range(model.size()):
                par = model[k]
                if par.is_free():
                
                    # Set parameter name
                    name = model_name+"_"+par.name()
                    
                    # Append value
                    colnames.append(name)
                    values[name] = par.value()
                    
                    # Append error
                    name = "Unc_"+name
                    colnames.append(name)
                    values[name] = par.error()
        
        # Bundle together results
        result = {'colnames': colnames, 'values': values}
        
        # Return
        return result


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates TS distribution.
    """
    # Create instance of application
    app = cstsdist(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
    