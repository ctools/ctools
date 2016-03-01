#! /usr/bin/env python
# ==========================================================================
# This script generates the TS distribution for a particular model based
# on Monte-Carlo simulations.
#
# Copyright (C) 2011-2016 Juergen Knoedlseder
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


# ============== #
# cstsdist class #
# ============== #
class cstsdist(ctools.cscript):
    """
    This class implements the TS distribution generation script. It derives
    from the ctools.cscript class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "cstsdist"
        self.version = "1.0.0"

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
            pars.append(gammalib.GApplicationPar("inobs","f","h","NONE","","","Input event list, counts cube, or observation definition XML file"))
            pars.append(gammalib.GApplicationPar("inmodel","f","a","$CTOOLS/share/models/crab.xml","","","Input model XML file"))
            pars.append(gammalib.GApplicationPar("srcname","s","a","Crab","","","Source name"))
            pars.append(gammalib.GApplicationPar("expcube","s","a","NONE","","","Input exposure cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("psfcube","s","a","NONE","","","Input PSF cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("bkgcube","s","a","NONE","","","Input background cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","prod2","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","South_0.5h","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("deadc","r","h","0.95","","","Deadtime correction factor"))
            pars.append(gammalib.GApplicationPar("edisp","b","h","no","","","Apply energy dispersion?"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","ts.dat","","","Output Test Statistics distribution file"))
            pars.append(gammalib.GApplicationPar("ntrials","i","a","10","","","Number of trials"))
            pars.append(gammalib.GApplicationPar("ra","r","a","83.63","0","360","RA of pointing (deg)"))
            pars.append(gammalib.GApplicationPar("dec","r","a","22.01","-90","90","Dec of pointing (deg)"))
            pars.append(gammalib.GApplicationPar("emin","r","a","0.1","","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","a","100.0","","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","0","","","Number of energy bins (0=unbinned)"))
            pars.append(gammalib.GApplicationPar("tmin","r","h","0.0","","","Start time (MET in s)"))
            pars.append(gammalib.GApplicationPar("tmax","r","a","1800.0","","","Observation duration (in s)"))
            pars.append(gammalib.GApplicationPar("npix","i","a","200","","","Number of pixels for binned"))
            pars.append(gammalib.GApplicationPar("binsz","r","a","0.05","","","Pixel size for binned (deg/pixel)"))
            pars.append(gammalib.GApplicationPar("rad","r","h","5.0","","","Radius of ROI (deg)"))         
            pars.append(gammalib.GApplicationPar("pattern","s","h","single","","","Observation pattern (single/four)"))
            pars.append(gammalib.GApplicationPar("offset","r","h","1.5","","","Observation pattern offset (deg)"))
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

            # Check for requested pattern and use above
            # observation parameters to set wobble pattern
            self.m_pattern = self["pattern"].string()
            if self.m_pattern == "four":
                self.obs = self.set_obs()

        # Get source name
        self.m_srcname = self["srcname"].string()

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

        # Set models if we have none
        if self.obs.models().size() == 0:
            self.obs.models(self["inmodel"].filename())

        # Get other parameters
        self.m_outfile = self["outfile"].filename()
        self.m_ntrials = self["ntrials"].integer()
        self.m_edisp   = self["edisp"].boolean()
        self.m_offset  = self["offset"].real()

        # Set some fixed parameters
        self.m_log   = False                    # Logging in client tools
        self.m_debug = self["debug"].boolean()  # Debugging in client tools

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
            result = self.trial(seed, full_model, bkg_model)

            # Write out result immediately
            if seed == 0:
                file   = open(self.m_outfile.url(), 'w')
                writer = csv.DictWriter(file, result['colnames'])
                headers = {}
                for n in result['colnames']:
                    headers[n] = n
                writer.writerow(headers)
                #writer.writerow(dict((_,_) for _ in result['colnames']))
            else:
                file = open(self.m_outfile.url(), 'a')
            writer = csv.DictWriter(file, result['colnames'])
            writer.writerow(result['values'])
            file.close()

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
            raise RuntimeError(msg)

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

    def set_obs(self):
        """
        Returns an observation container with a set of CTA observations.

        Keywords:
        """

        # Setup observation definition list
        obsdeflist = obsutils.set_obs_patterns(self.m_pattern,
                                               ra=self["ra"].real(),
                                               dec=self["dec"].real(),
                                               offset=self["offset"].real())

        # Create list of observations
        obs = obsutils.set_obs_list(obsdeflist,
                                    tstart=self["tmin"].real(),
                                    duration=self["tmax"].real()-self["tmin"].real(),
                                    deadc=self["deadc"].real(),
                                    emin=self["emin"].real(),
                                    emax=self["emax"].real(),
                                    rad=self["rad"].real(),
                                    irf=self["irf"].string(), caldb=self["caldb"].string())

        # Return observation container
        return obs

    def trial(self, seed, full_model, bkg_model):
        """
        Create the TS for a single trial.

        Parameters:
         seed - Random number generator seed
        """
        # Write header
        if self.logExplicit():
            self.log.header2("Trial "+str(seed+1))

        # Simulate events
        sim = obsutils.sim(self.obs,
                           nbins=self.m_enumbins,
                           seed=seed,
                           binsz=self.m_binsz,
                           npix=self.m_npix,
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
        like_bgm   = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
        result_bgm = like_bgm.obs().models().copy()
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
        sim.models(full_model)
        like_all   = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
        result_all = like_all.obs().models().copy()
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
            self.log(result_all[self.m_srcname]["Prefactor"].value())
            self.log("+/-")
            self.log(result_all[self.m_srcname]["Prefactor"].error())
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
