#! /usr/bin/env python
# ==========================================================================
# This script computes the CTA sensitivity.
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
#from ctools import obsutils
from cscripts import obsutils
import sys
import csv
import math


# ============ #
# cssens class #
# ============ #
class cssens(ctools.cscript):
    """
    This class computes the CTA sensitivity for a number of energy bins using
    ctlike. Spectra are fitted in narrow energy bins to simulated data,
    and the flux level is determined that leads to a particular significance.

    The significance is determined using the Test Statistic, defined as twice
    the likelihood difference between fitting with and  without the test source.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "cssens"
        self.version = "1.0.0"

        # Initialise some parameters
        self.obs    = None
        self.m_ra   = None
        self.m_dec  = None
        
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
        self.logFileOpen()
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
            pars.append(gammalib.GApplicationPar("srcname","s","a","Crab","","","Source name"))
            pars.append(gammalib.GApplicationPar("offset","r","h","0.0","","","Source offset angle (deg)"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","prod2","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","South_50h","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("deadc","r","h","0.95","0","1","Deadtime correction factor"))
            pars.append(gammalib.GApplicationPar("outfile","f","h","sensitivity.dat","","","Output file name"))
            pars.append(gammalib.GApplicationPar("duration","r","a","180000.0","","","Effective exposure time (s)"))
            pars.append(gammalib.GApplicationPar("rad","r","a","5.0","","","Radius of ROI (deg)"))
            pars.append(gammalib.GApplicationPar("emin","r","a","0.020","","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","a","200.0","","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("bins","i","a","21","","","Number of energy bins for differential sensitivity computation"))
            pars.append(gammalib.GApplicationPar("enumbins","i","h","0","","","Number of energy bins (0=unbinned)"))
            pars.append(gammalib.GApplicationPar("npix","i","h","200","","","Number of pixels for binned"))
            pars.append(gammalib.GApplicationPar("binsz","r","h","0.05","","","Pixel size for binned (deg/pixel)"))
            pars.append(gammalib.GApplicationPar("type","s","h","Differential","Differential|Integral","","Sensitivity type"))
            pars.append(gammalib.GApplicationPar("sigma","r","h","5.0","","","Significance threshold"))
            pars.append(gammalib.GApplicationPar("max_iter","i","h","50","","","Maximum number of iterations"))
            pars.append(gammalib.GApplicationPar("num_avg","i","h","3","","","Number of iterations for sliding average"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","cssens.log","","","Log filename"))
            pars.save(parfile)
        
        # Return
        return
        
    def get_parameters(self):
        """
        Get user parameters from parfile.
        """
        # Set observation if not done before
        if self.obs == None or self.obs.size() == 0:
            self.obs = self.set_obs()

        # Set models if we have none
        if self.obs.models().size() == 0:
            self.obs.models(self["inmodel"].filename())

        # Get source name
        self.m_srcname = self["srcname"].string()

        # Read further parameters
        self.m_outfile = self["outfile"].filename()
        self.m_emin    = self["emin"].real()
        self.m_emax    = self["emax"].real()
        self.m_bins    = self["bins"].integer()

        # Read parameters for binned if requested
        self.m_enumbins = self["enumbins"].integer()
        if not self.m_enumbins == 0:
            self.m_npix  = self["npix"].integer()
            self.m_binsz = self["binsz"].real()
        else:
            # Set dummy values (required by obsutils)
            self.m_npix  = 0
            self.m_binsz = 0.0

        # Read remaining parameters
        self.m_ts_thres = self["sigma"].real()*self["sigma"].real()
        self.m_max_iter = self["max_iter"].integer()
        self.m_num_avg  = self["num_avg"].integer()
        self.m_type     = self["type"].string()
        
        # Set some fixed parameters
        self.m_log   = False                   # Logging in client tools
        self.m_debug = self["debug"].boolean() # Debugging in client tools

        # Derive some parameters
        self.m_ebounds = gammalib.GEbounds(self.m_bins, \
                                           gammalib.GEnergy(self.m_emin, "TeV"), \
                                           gammalib.GEnergy(self.m_emax, "TeV"))
        
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
        
        # Initialise script
        colnames = ['loge', 'emin', 'emax', 'crab_flux', 'photon_flux', \
                    'energy_flux', 'sensitivity']
        results  = []

        # Initialise models
        full_model, bkg_model = self.set_models()

        # Write models into logger
        if self.logTerse():
            self.log("\n")
            self.log.header1("Models")
            self.log.header2("Background model")
            self.log(str(bkg_model))
            self.log("\n\n")
            self.log.header2("Full model")
            self.log(str(full_model))
            self.log("\n")

        # Write heder
        if self.logTerse():
            self.log("\n")
            self.log.header1("Sensitivity determination")
            self.log.parformat("Type")
            self.log(self.m_type)
            self.log("\n")

        # Loop over energy bins
        for ieng in range(self.m_ebounds.size()):
        
            # Set energies
            if self.m_type == "Differential":
                emin  = self.m_ebounds.emin(ieng)
                emax  = self.m_ebounds.emax(ieng)
            elif self.m_type == "Integral":
                emin  = self.m_ebounds.emin(ieng)
                emax  = self.m_ebounds.emax()
            else:
                msg = "Invalid sensitivity type \""+self.m_type+"\" encountered."+ \
                      " Either use \"Differential\" or \"Integral\"."
                raise gammalib.GException.invalid_value("cssens", msg)
            
            # Determine sensitivity
            result = self.get_sensitivity(self.obs, emin, emax, \
                                          bkg_model, full_model)
            
            # Write results
            if ieng == 0:
                f      = open(self.m_outfile, 'w')
                writer = csv.DictWriter(f, colnames)
                writer.writerow(dict((_,_) for _ in colnames))
                writer.writerow(result)
                f.close()
            else:
                f = open(self.m_outfile, 'a')
                writer = csv.DictWriter(f, colnames)
                writer.writerow(result)
                f.close()
            
            # Append results
            results.append(result)
        
        # Return
        return

    def set_obs(self, lpnt=0.0, bpnt=0.0, emin=0.1, emax=100.0):
        """
        Returns an observation container.
        
        Keywords:
         lpnt - Galactic longitude of pointing [deg] (default: 0.0)
         bpnt - Galactic latitude of pointing [deg] (default: 0.0)
         emin - Minimum energy [TeV] (default: 0.1)
         emax - Maximum energy [TeV] (default: 100.0)
        """
        # If an observation was provided on input then load it from XML file
        filename = self["inobs"].filename()
        if filename != "NONE" and filename != "":
            obs = self.get_observations()

        # ... otherwise allocate a single observation
        else:

            # Read relevant user parameters
            caldb    = self["caldb"].string()
            irf      = self["irf"].string()
            deadc    = self["deadc"].real()
            duration = self["duration"].real()
            rad      = self["rad"].real()

            # Allocate observation container
            obs = gammalib.GObservations()
    
            # Set single pointing
            pntdir = gammalib.GSkyDir()
            pntdir.lb_deg(lpnt, bpnt)
        
            # Create CTA observation
            run = obsutils.set_obs(pntdir, caldb=caldb, irf=irf, \
                                   duration=duration, deadc=deadc, \
                                   emin=emin, emax=emax, rad=rad)
        
            # Append observation to container
            obs.append(run)

            # Set source position
            offset     = self["offset"].real()
            pntdir.lb_deg(lpnt, bpnt+offset)
            self.m_ra  = pntdir.ra_deg()
            self.m_dec = pntdir.dec_deg()
    
        # Return observation container
        return obs

    def set_obs_ebounds(self, emin, emax):
        """
        Set energy boundaries for observation in container.
        
        Parameters:
         emin - Minimum energy
         emax - Maximum energy
        """
        # Loop over all observations in container
        for obs in self.obs:
        
            # Set energy boundaries
            ebounds = gammalib.GEbounds(emin, emax)
            obs.events().ebounds(ebounds)
        
        # Return
        return

    def set_models(self, fitpos=False, fitspec=False):
        """
        Set full and background model.
        """
        # Retrieve full model from observation container
        full_model = self.obs.models().copy()

        # Get source model
        model = full_model[self.m_srcname]
        
        # Check that model has a Prefactor
        if not model.has_par("Prefactor"):
            msg = "Model \""+self.m_srcname+"\" has no parameter \"Prefactor\"."+ \
                  " Only spectral models with a \"Prefactor\" parameter are supported."
            raise gammalib.GException.invalid_value("cssens", msg)

        # Set source position
        if self.m_ra != None and self.m_dec != None:
            if model.has_par("RA") and model.has_par("DEC"):
                model["RA"].value(self.m_ra)
                model["DEC"].value(self.m_dec)

        # Fit or fix spatial parameters
        if fitpos:
            if model.has_par("RA"):
                model["RA"].free()
            if model.has_par("DEC"):
                model["DEC"].free()
            if model.has_par("Sigma"):
                model["Sigma"].free()
            if model.has_par("Radius"):
                model["Radius"].free()
            if model.has_par("Width"):
                model["Width"].free()
        else:
            if model.has_par("RA"):
                model["RA"].fix()
            if model.has_par("DEC"):
                model["DEC"].fix()
            if model.has_par("Sigma"):
                model["Sigma"].fix()
            if model.has_par("Radius"):
                model["Radius"].fix()
            if model.has_par("Width"):
                model["Width"].fix()

        # Fit or fix spectral parameters
        if fitspec:
            if model.has_par("Index"):
                model["Index"].free()
            if model.has_par("Cutoff"):
                model["Cutoff"].free()
        else:
            if model.has_par("Index"):
                model["Index"].fix()
            if model.has_par("Cutoff"):
                model["Cutoff"].fix()

        # Create background model
        bkg_model = full_model.copy()
        bkg_model.remove(self.m_srcname)

        # Return models
        return full_model, bkg_model

    def get_sensitivity(self, obs, emin, emax, bkg_model, full_model):
        """
        Determine sensitivity for given observations.
        
        Parameters:
         obs        - Observation container
         emin       - Minimum energy for fitting and flux computation
         emax       - Maximum energy for fitting and flux computation
         bkg_model  - Background model
         full_model - Source model
        """
        # Set TeV->erg conversion factor
        tev2erg = 1.6021764

        # Set energy boundaries
        self.set_obs_ebounds(emin, emax)
        
        # Determine energy boundaries from first observation in the container
        loge      = math.log10(math.sqrt(emin.TeV()*emax.TeV()))
        e_mean    = math.pow(10.0, loge)
        erg_mean  = e_mean * tev2erg

        # Compute Crab unit (this is the factor with which the Prefactor needs
        # to be multiplied to get 1 Crab
        crab_flux = self.get_crab_flux(emin, emax)
        src_flux  = full_model[self.m_srcname].spectral().flux(emin, emax)
        crab_unit = crab_flux/src_flux

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header2("Energies: "+str(emin)+" - "+str(emax))
            self.log.parformat("Crab flux")
            self.log(crab_flux)
            self.log(" ph/cm2/s\n")
            self.log.parformat("Source model flux")
            self.log(src_flux)
            self.log(" ph/cm2/s\n")
            self.log.parformat("Crab unit factor")
            self.log(crab_unit)
            self.log("\n")

        # Initialise loop
        crab_flux_value   = []
        photon_flux_value = []
        energy_flux_value = []
        sensitivity_value = []
        iter              = 0
        test_crab_flux    = 0.1 # Initial test flux in Crab units (100 mCrab)

        # Loop until we break
        while True:

            # Update iteration counter
            iter += 1

            # Write header
            if self.logExplicit():
                self.log.header2("Iteration "+str(iter))

            # Set source model. crab_prefactor is the Prefactor that
            # corresponds to 1 Crab
            src_model      = full_model.copy()
            crab_prefactor = src_model[self.m_srcname]['Prefactor'].value() * crab_unit
            src_model[self.m_srcname]['Prefactor'].value(crab_prefactor * test_crab_flux)
            obs.models(src_model)

            # Simulate events
            sim = obsutils.sim(obs, nbins=self.m_enumbins, seed=iter, \
                               binsz=self.m_binsz, npix=self.m_npix, \
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
            sim.models(bkg_model)
            like       = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
            result_bgm = like.obs().models().copy()
            LogL_bgm   = like.opt().value()
            npred_bgm  = like.obs().npred()

            # Assess quality based on a comparison between Npred and Nevents
            quality_bgm = npred_bgm-nevents

            # Write background fit results
            if self.logExplicit():
                self.log.header3("Background model fit")
                self.log.parformat("log likelihood")
                self.log(LogL_bgm)
                self.log("\n")
                self.log.parformat("Number of predicted events")
                self.log(npred_bgm)
                self.log("\n")
                self.log.parformat("Fit quality")
                self.log(quality_bgm)
                self.log("\n")

            # Start over if the fit quality was bad
            if abs(quality_bgm) > 3.0:
                if self.logExplicit():
                    self.log("Fit quality outside required range. Start over.\n")
                continue

            # Write model fit results
            if self.logExplicit():
                for model in result_bgm:
                    self.log.parformat("Model")
                    self.log(model.name())
                    self.log("\n")
                    for par in model:
                        self.log(str(par)+"\n")
            
            # Fit background and test source
            sim.models(src_model)
            like       = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
            result_all = like.obs().models().copy()
            LogL_all   = like.opt().value()
            npred_all  = like.obs().npred()
            ts         = 2.0*(LogL_bgm-LogL_all)

            # Assess quality based on a comparison between Npred and Nevents
            quality_all = npred_all-nevents

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
                self.log.parformat("Fit quality")
                self.log(quality_all)
                self.log("\n")
                #
                for model in result_all:
                    self.log.parformat("Model")
                    self.log(model.name())
                    self.log("\n")
                    for par in model:
                        self.log(str(par)+"\n")

            # Start over if the fit quality was bad
            if abs(quality_all) > 3.0:
                if self.logExplicit():
                    self.log("Fit quality outside required range. Start over.\n")
                continue

            # Start over if TS was non-positive
            if ts <= 0.0:
                if self.logExplicit():
                    self.log("Non positive TS. Start over.\n")
                continue

            # Get fitted Crab, photon and energy fluxes
            crab_flux     = result_all[self.m_srcname]['Prefactor'].value() / crab_prefactor
            crab_flux_err = result_all[self.m_srcname]['Prefactor'].error() / crab_prefactor
            photon_flux   = result_all[self.m_srcname].spectral().flux(emin, emax)
            energy_flux   = result_all[self.m_srcname].spectral().eflux(emin, emax)

            # Compute differential sensitivity in unit erg/cm2/s
            energy      = gammalib.GEnergy(e_mean, "TeV")
            time        = gammalib.GTime()
            sensitivity = result_all[self.m_srcname].spectral().eval(energy, time) * \
                          erg_mean*erg_mean * 1.0e6

            # Compute flux correction factor based on average TS
            correct = 1.0
            if ts > 0:
                correct = math.sqrt(self.m_ts_thres/ts)
            
            # Compute extrapolated fluxes
            crab_flux   = correct * crab_flux
            photon_flux = correct * photon_flux
            energy_flux = correct * energy_flux
            sensitivity = correct * sensitivity
            crab_flux_value.append(crab_flux)
            photon_flux_value.append(photon_flux)
            energy_flux_value.append(energy_flux)
            sensitivity_value.append(sensitivity)
            
            # Write background and test source fit results
            if self.logExplicit():
                self.log.parformat("Photon flux")
                self.log(photon_flux)
                self.log(" ph/cm2/s\n")
                self.log.parformat("Energy flux")
                self.log(energy_flux)
                self.log(" erg/cm2/s\n")
                self.log.parformat("Crab flux")
                self.log(crab_flux*1000.0)
                self.log(" mCrab\n")
                self.log.parformat("Differential sensitivity")
                self.log(sensitivity)
                self.log(" erg/cm2/s\n")
                for model in result_all:
                    self.log.parformat("Model")
                    self.log(model.name())
                    self.log("\n")
                    for par in model:
                        self.log(str(par)+"\n")
            elif self.logTerse():
                self.log.parformat("Iteration "+str(iter))
                self.log("TS=")
                self.log(ts)
                self.log(" ")
                self.log("corr=")
                self.log(correct)
                self.log("  ")
                self.log(photon_flux)
                self.log(" ph/cm2/s = ")
                self.log(energy_flux)
                self.log(" erg/cm2/s = ")
                self.log(crab_flux*1000.0)
                self.log(" mCrab = ")
                self.log(sensitivity)
                self.log(" erg/cm2/s\n")
            
            # Compute sliding average of extrapolated fitted prefactor,
            # photon and energy flux. This damps out fluctuations and
            # improves convergence
            crab_flux   = 0.0
            photon_flux = 0.0
            energy_flux = 0.0
            sensitivity = 0.0
            num         = 0.0
            for k in range(self.m_num_avg):
                inx = len(crab_flux_value) - k - 1
                if inx >= 0:
                    crab_flux   += crab_flux_value[inx]
                    photon_flux += photon_flux_value[inx]
                    energy_flux += energy_flux_value[inx]
                    sensitivity += sensitivity_value[inx]
                    num      += 1.0
            crab_flux   /= num
            photon_flux /= num
            energy_flux /= num
            sensitivity /= num
            
            # Compare average flux to last average
            if iter > self.m_num_avg:
                if test_crab_flux > 0:
                    ratio = crab_flux/test_crab_flux
                    
                    # We have 2 convergence criteria:
                    # 1. The average flux does not change
                    # 2. The flux correction factor is small
                    if ratio   >= 0.99 and ratio   <= 1.01 and \
                       correct >= 0.9  and correct <= 1.1:
                        if self.logTerse():
                            self.log(" Converged ("+str(ratio)+")\n")
                        break
                else:
                    if self.logTerse():
                        self.log(" Flux is zero.\n")
                    break
            
            # Use average for next iteration
            test_crab_flux = crab_flux
            
            # Exit loop if number of trials exhausted
            if (iter >= self.m_max_iter):
                if self.logTerse():
                    self.log(" Test ended after "+str(self.m_max_iter)+" iterations.\n")
                break

        # Write fit results
        if self.logTerse():
            self.log.header3("Fit results")
            self.log.parformat("Test statistics")
            self.log(ts)
            self.log("\n")
            self.log.parformat("Photon flux")
            self.log(photon_flux)
            self.log(" ph/cm2/s\n")
            self.log.parformat("Energy flux")
            self.log(energy_flux)
            self.log(" erg/cm2/s\n")
            self.log.parformat("Crab flux")
            self.log(crab_flux*1000.0)
            self.log(" mCrab\n")
            self.log.parformat("Differential sensitivity")
            self.log(sensitivity)
            self.log(" erg/cm2/s\n")
            self.log.parformat("Number of simulated events")
            self.log(nevents)
            self.log("\n")
            self.log.header3("Background and test source model fitting")
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
            
        # Store result
        result = {'loge': loge, 'emin': emin.TeV(), 'emax': emax.TeV(), \
                  'crab_flux': crab_flux, 'photon_flux': photon_flux, \
                  'energy_flux': energy_flux, \
                  'sensitivity': sensitivity}
        
        # Return result
        return result

    def get_crab_flux(self, emin, emax):
        """
        Determine the Crab photon flux in a given energy interval.
        """
        # Set Crab spectral model
        crab = gammalib.GModelSpectralPlaw(5.7e-16, -2.48, gammalib.GEnergy(0.3, "TeV"))

        # Determine flux
        flux = crab.flux(emin, emax)

        # Return flux
        return flux


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    CTA ctlike-based sensitivity estimator.
    """
    # Create instance of application
    app = cssens(sys.argv)
    
    # Run application
    app.execute()
    