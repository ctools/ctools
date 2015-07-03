#! /usr/bin/env python
# ==========================================================================
# This script generates the pull distribution for all model parameters.
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
from cscripts import obsutils
import sys
import csv


# ============ #
# cspull class #
# ============ #
class cspull(ctools.cscript):
    """
    This class implements the pull distribution generation script. It derives
    from the ctools.cscript class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "cspull"
        self.version = "1.0.0"
        
        # Initialise some members
        self.obs           = None
        self.model         = None
        self.m_inmodel     = None
        self.m_edisp       = False
        self.m_profile     = False
        self.m_exposure    = None
        self.m_psfcube     = None
        self.m_bckcube     = None
        self.m_stackmodels = None
        
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
            print("Parfile "+parfile+" not found. Create default parfile.")
            
            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("inobs","f","h","NONE","","","Event list, counts cube, or observation definition file"))
            pars.append(gammalib.GApplicationPar("inmodel","f","a","$CTOOLS/share/models/crab.xml","","","Source model"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","pull.dat","","","Output file name"))
            #pars.append(gammalib.GApplicationPar("expcube","s","a","NONE","","","Exposure cube file (only needed for stacked analysis)"))
            #pars.append(gammalib.GApplicationPar("psfcube","s","a","NONE","","","PSF cube file (only needed for stacked analysis)"))
            #pars.append(gammalib.GApplicationPar("bkgcube","s","a","NONE","","","Background cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","prod2","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","South_50h","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("deadc","r","h","0.95","0","1","Deadtime correction factor"))
            pars.append(gammalib.GApplicationPar("edisp","b","h","no","","","Apply energy dispersion?"))
            pars.append(gammalib.GApplicationPar("profile","b","h","no","","","Use likelihood profile method for errors?"))
            pars.append(gammalib.GApplicationPar("ntrials","i","a","10","","","Number of trials"))
            pars.append(gammalib.GApplicationPar("ra","r","a","83.6331","0","360","RA of pointing (deg)"))
            pars.append(gammalib.GApplicationPar("dec","r","a","22.0145","-90","90","Dec of pointing (deg)"))
            pars.append(gammalib.GApplicationPar("emin","r","a","0.1","","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","a","100.0","","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","0","","","Number of energy bins (0=unbinned)"))
            pars.append(gammalib.GApplicationPar("tmin","r","h","0.0","","","Start time (MET in s)"))
            pars.append(gammalib.GApplicationPar("tmax","r","a","1800.0","","","Duration (in s)"))
            pars.append(gammalib.GApplicationPar("npix","i","a","200","","","Number of pixels for binned"))
            pars.append(gammalib.GApplicationPar("binsz","r","a","0.05","","","Pixel size for binned (deg/pixel)"))
            pars.append(gammalib.GApplicationPar("rad","r","h","5.0","0","180","Radius of ROI (deg)"))         
            pars.append(gammalib.GApplicationPar("pattern","s","h","single","","","Observation pattern (single/four)"))
            pars.append(gammalib.GApplicationPar("offset","r","h","1.5","","","Observation pattern offset (deg)"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","cspull.log","","","Log filename"))
            pars.save(parfile)
        
        # Return
        return
        
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters
        
        # Set observation if not done before
        if self.obs == None or self.obs.size() == 0:
            self.obs = self.get_observations()
            
            # Check for requested pattern and use above
            # observation parameters to set wobble pattern
            self.m_pattern = self["pattern"].string()
            if self.m_pattern == "four":
                self.obs = self.set_obs()
            
        # Get number of energy bins
        self.m_enumbins = self["enumbins"].integer()
   
        # Read parameters for binned if requested
        if not self.m_enumbins == 0:
            self.m_npix  = self["npix"].integer()
            self.m_binsz = self["binsz"].real()
        else:
            # Set dummy values (required by obsutils)
            self.m_npix  = 0
            self.m_binsz = 0.0
            
        # Set models if we have none
        if self.obs.models().size() == 0:
            self.obs.models(self["inmodel"].filename())
         
        # Read other parameters    
        self.m_outfile = self["outfile"].filename()
        self.m_ntrials = self["ntrials"].integer()   
        self.m_edisp   = self["edisp"].boolean()
        self.m_offset  = self["offset"].real()   
        self.m_profile = self["profile"].boolean()

        # Set some fixed parameters
        self.m_log     = False                     # Logging in client tools
        self.m_chatter = self["chatter"].integer() # Chatter level 
        self.m_debug   = self["debug"].boolean()   # Debugging in client tools

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

        
    def set_stacked_irf(self):
        """
        For stacked analysis prepare stacked irfs
        """
        # Write header into logger
        if self.logTerse():
            self.log("\n")
            self.log.header1("Compute stacked response")

        # Get stacked exposure
        expcube = ctools.ctexpcube(self.obs)
        expcube["incube"]   = "NONE"
        expcube["usepnt"]   = True
        expcube["ebinalg"]  = "LOG"
        expcube["binsz"]    = self.m_binsz
        expcube["nxpix"]    = self.m_npix
        expcube["nypix"]    = self.m_npix
        expcube["enumbins"] = self.m_enumbins
        expcube["emin"]     = self["emin"].real()
        expcube["emax"]     = self["emax"].real()
        expcube["coordsys"] = "CEL"
        expcube["proj"]     = "TAN" 
        expcube.run()

        # Notify exposure computation
        if self.logTerse():
            self.log("Computed exposure cube\n")
        
        # Get stacked Psf
        psfcube = ctools.ctpsfcube(self.obs)
        psfcube["incube"]   = "NONE"
        psfcube["usepnt"]   = True
        psfcube["ebinalg"]  = "LOG"
        psfcube["binsz"]    = self.m_binsz*10.0
        psfcube["nxpix"]    = self.m_npix/10
        psfcube["nypix"]    = self.m_npix/10
        psfcube["enumbins"] = self.m_enumbins
        psfcube["emin"]     = self["emin"].real()
        psfcube["emax"]     = self["emax"].real()
        psfcube["coordsys"] = "CEL"
        psfcube["proj"]     = "TAN"      
        psfcube.run()

        # Notify Psf computation
        if self.logTerse():
            self.log("Computed Psf cube\n")
        
        # Get stacked background
        bkgcube = ctools.ctbkgcube(self.obs)
        bkgcube["incube"]   = "NONE"
        bkgcube["usepnt"]   = True
        bkgcube["ebinalg"]  = "LOG"
        bkgcube["binsz"]    = self.m_binsz
        bkgcube["nxpix"]    = self.m_npix
        bkgcube["nypix"]    = self.m_npix
        bkgcube["enumbins"] = self.m_enumbins
        bkgcube["emin"]     = self["emin"].real()
        bkgcube["emax"]     = self["emax"].real()
        bkgcube["coordsys"] = "CEL"
        bkgcube["proj"]     = "TAN"
        bkgcube.run()

        # Notify background cube computation
        if self.logTerse():
            self.log("Computed background cube\n")
        
        # Store results
        self.m_exposure    = expcube.expcube().copy()
        self.m_psfcube     = psfcube.psfcube().copy()
        self.m_bckcube     = bkgcube.bkgcube().copy()
        self.m_stackmodels = bkgcube.models().copy()
        
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
            self.log.header1("Observation(s)")
            self.log(str(self.obs))
            self.log("\n")

        # If several observations and binned: prepare stacked irfs
        if self.obs.size() > 1 and self.m_enumbins > 0:
            self.set_stacked_irf()

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Generate pull distribution")

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
    

    def set_obs(self):
        """
        Returns an observation container with a set of CTA observations.
        
        Keywords:
        """
        
        # Setup observation definition list
        obsdeflist = obsutils.set_obs_patterns(self.m_pattern, \
                                               ra=self["ra"].real(), dec=self["dec"].real(), \
                                               offset=self["offset"].real())
        
        # Create list of observations
        obs = obsutils.set_obs_list(obsdeflist, \
                                    tstart=self["tmin"].real(), duration=self["tmax"].real()-self["tmin"].real(), \
                                    deadc=self["deadc"].real(), \
                                    emin=self["emin"].real(), emax=self["emax"].real(), \
                                    rad=self["rad"].real(), \
                                    irf=self["irf"].string(), caldb=self["caldb"].string())
    
        # Return observation container
        return obs


    def trial(self, seed):
        """
        Create the pull for a single trial.
        
        Parameters:
         seed - Random number generator seed
        """
        
        # Write header
        if self.logNormal():
            self.log.header2("Trial "+str(seed+1))

        # Simulate events
        obs = obsutils.sim(self.obs, \
                           nbins=self.m_enumbins, \
                           seed=seed, \
                           binsz=self.m_binsz, \
                           npix=self.m_npix, \
                           edisp=self.m_edisp, \
                           log=self.m_log, \
                           debug=self.m_debug, \
                           chatter=self.m_chatter)
        
        # If stacked, add stacked responses and model
        if self.obs.size() > 0 and self.m_enumbins > 0:
            obs[0].response(self.m_exposure, self.m_psfcube, self.m_bckcube)
            obs.models(self.m_stackmodels)
        
        # Determine number of events in simulation
        nevents = 0.0
        for run in obs:
            nevents += run.events().number()

        # Write simulation results
        if self.logNormal():
            self.log.header3("Simulation")
            self.log.parformat("Number of simulated events")
            self.log(nevents)
            self.log("\n")

        # Fit model
        if self.m_profile:
            models = obs.models()
            for i in range(models.size()):
                model_name = models[i].name()
                like       = obsutils.cterror(obs, model_name, \
                                              log=self.m_log, \
                                              debug=self.m_debug,
                                              chatter=self.m_chatter)
        else:
            like = obsutils.fit(obs, edisp=self.m_edisp, \
                                log=self.m_log, \
                                debug=self.m_debug, \
                                chatter=self.m_chatter)

        # Store results
        logL   = like.opt().value()
        npred  = like.obs().npred()
        models = like.obs().models()

        # Write result header
        if self.logNormal():
            self.log.header3("Pulls")
        
        # Gather results
        colnames = []
        values   = {}
        colnames.append("LogL")
        colnames.append("Sim_Events")
        colnames.append("Npred_Events")
        values["LogL"]         = logL
        values["Sim_Events"]   = nevents
        values["Npred_Events"] = npred
        for i in range(models.size()):
            model      = models[i]
            model_name = model.name()
            for k in range(model.size()):
                par = model[k]
                if par.is_free():
                
                    # Set parameter name
                    name = model_name+"_"+par.name()
                    
                    # Append parameter, Pull_parameter and Unc_parameter
                    colnames.append(name)
                    colnames.append("Pull_"+name)
                    colnames.append("Unc_"+name)
                
                    # Compute pull
                    fitted_value = par.value()
                    real_value   = self.obs.models()[i][k].value()
                    error        = par.error()
                    if error != 0.0:
                        pull = (fitted_value - real_value) / error
                    else:
                        pull = 99.0
                        
                    # Store results
                    values[name] = fitted_value
                    values["Pull_"+name] = pull
                    values["Unc_"+name] = error

                    # Write result
                    if self.logNormal():
                        self.log.parformat(name)
                        self.log(pull)
                        self.log(" (")
                        self.log(fitted_value)
                        self.log(" +/- ")
                        self.log(error)
                        self.log(")\n")
        
        # Bundle together results
        result = {'colnames': colnames, 'values': values}
        
        # Return
        return result


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates pull distribution for a source model.
    """
    # Create instance of application
    app = cspull(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
    
