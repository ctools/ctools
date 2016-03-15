#! /usr/bin/env python
# ==========================================================================
# Generates the TS distribution for a particular model
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
    Generates Test Statistic distribution for a model.
    
    
    This class implements the TS distribution generation script. It derives
    from the ctools.cscript class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = "cstsdist"
        self._version = "1.1.0"

        # Initialise some members
        self._obs         = gammalib.GObservations()
        self._pattern     = "single"
        self._srcname     = ""
        self._enumbins    = 0
        self._npix        = 0
        self._binsz       = 0.0
        self._outfile     = gammalib.GFilename("ts.dat")
        self._ntrials     = 10
        self._edisp       = False
        self._debug       = False
        self._log_clients = False        

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

        # Set observation if not done before
        if self._obs == None or self._obs.size() == 0:
            self._obs = self._get_observations()

            # Check for requested pattern and use above
            # observation parameters to set wobble pattern
            self._pattern = self["pattern"].string()
            if self._pattern == "four":
                self._obs = self._set_obs()

        # Get source name
        self._srcname = self["srcname"].string()

        # Get number of energy bins
        self._enumbins = self["enumbins"].integer()

        # Read parameters for binned if requested
        if not self._enumbins == 0:
            self._npix  = self["npix"].integer()
            self._binsz = self["binsz"].real()

        # Set models if we have none
        if self._obs.models().size() == 0:
            self._obs.models(self["inmodel"].filename())

        # Get other parameters
        self._outfile = self["outfile"].filename()
        self._ntrials = self["ntrials"].integer()
        self._edisp   = self["edisp"].boolean()
        self._debug   = self["debug"].boolean()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Return
        return

    def _set_obs(self):
        """
        Set an observation container.

        Returns:
            Observation container.
        """

        # Setup observation definition list
        obsdeflist = obsutils.set_obs_patterns(self._pattern,
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
                                    irf=self["irf"].string(),
                                    caldb=self["caldb"].string())

        # Return observation container
        return obs

    def _set_models(self, fitspat=False, fitspec=False):
        """
        Set full model and background model.

        Kwargs:
            fitspat: Fit spatial parameter (default: False).
            fitspec: Fit spectral parameters (default: False).

        Returns:
            Tuple containing full model and background model.
        """
        # Retrieve full model from observation container
        full_model = self._obs.models().copy()

        # Get source model
        model = full_model[self._srcname]

        # Check that model has a Prefactor
        if not model.has_par("Prefactor"):
            msg = "Model \""+self._srcname+"\" has no parameter "+\
                  "\"Prefactor\". Only spectral models with a "+\
                  "\"Prefactor\" parameter are supported."
            raise RuntimeError(msg)

        # Fit or fix spatial parameters
        if fitspat:
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
        bkg_model.remove(self._srcname)

        # Return models
        return full_model, bkg_model

    def _trial(self, seed, full_model, bkg_model):
        """
        Create the TS for a single trial.

        Args:
            seed:       Random number generator seed
            full_model: Full model
            bkg_model:  Background model

        Returns:
            Result dictionary
        """
        # Write header
        if self._logExplicit():
            self._log.header2("Trial "+str(seed+1))

        # Simulate events
        sim = obsutils.sim(self._obs,
                           nbins=self._enumbins,
                           seed=seed,
                           binsz=self._binsz,
                           npix=self._npix,
                           log=self._log_clients,
                           debug=self._debug)

        # Determine number of events in simulation
        nevents = 0.0
        for run in sim:
            nevents += run.events().number()

        # Write simulation results
        if self._logExplicit():
            self._log.header3("Simulation")
            self._log.parformat("Number of simulated events")
            self._log(nevents)
            self._log("\n")

        # Fit background only
        sim.models(bkg_model)
        like_bgm   = obsutils.fit(sim,
                                  log=self._log_clients,
                                  debug=self._debug)
        result_bgm = like_bgm.obs().models().copy()
        LogL_bgm   = like_bgm.opt().value()
        npred_bgm  = like_bgm.obs().npred()

        # Write background fit results
        if self._logExplicit():
            self._log.header3("Background model fit")
            self._log.parformat("log likelihood")
            self._log(LogL_bgm)
            self._log("\n")
            self._log.parformat("Number of predicted events")
            self._log(npred_bgm)
            self._log("\n")
            for model in result_bgm:
                self._log.parformat("Model")
                self._log(model.name())
                self._log("\n")
                for par in model:
                    self._log(str(par)+"\n")

        # Fit background and test source
        sim.models(full_model)
        like_all   = obsutils.fit(sim,
                                  log=self._log_clients,
                                  debug=self._debug)
        result_all = like_all.obs().models().copy()
        LogL_all   = like_all.opt().value()
        npred_all  = like_all.obs().npred()
        ts         = 2.0*(LogL_bgm-LogL_all)

        # Write background and test source fit results
        if self._logExplicit():
            self._log.header3("Background and test source model fit")
            self._log.parformat("Test statistics")
            self._log(ts)
            self._log("\n")
            self._log.parformat("log likelihood")
            self._log(LogL_all)
            self._log("\n")
            self._log.parformat("Number of predicted events")
            self._log(npred_all)
            self._log("\n")
            for model in result_all:
                self._log.parformat("Model")
                self._log(model.name())
                self._log("\n")
                for par in model:
                    self._log(str(par)+"\n")

        # Write result
        elif self._logTerse():
            self._log.parformat("Trial "+str(seed))
            self._log("TS=")
            self._log(ts)
            self._log("  Prefactor=")
            self._log(result_all[self._srcname]["Prefactor"].value())
            self._log("+/-")
            self._log(result_all[self._srcname]["Prefactor"].error())
            self._log("\n")

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

        # Initialise models
        full_model, bkg_model = self._set_models()

        # Write models into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Models")
            self._log.header2("Background model")
            self._log(str(bkg_model))
            self._log("\n\n")
            self._log.header2("Full model")
            self._log(str(full_model))
            self._log("\n")

        # Write observation into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Observation")
            self._log(str(self._obs))
            self._log("\n")

        # Write models into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Test model")
            self._log(str(full_model))
            self._log("\n")

        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Generate TS distribution")

        # Loop over trials
        for seed in range(self._ntrials):

            # Make a trial
            result = self._trial(seed, full_model, bkg_model)

            # Write out result immediately
            if seed == 0:
                f      = open(self._outfile.url(), 'w')
                writer = csv.DictWriter(f, result['colnames'])
                headers = {}
                for n in result['colnames']:
                    headers[n] = n
                writer.writerow(headers)
            else:
                f = open(self._outfile.url(), 'a')
            writer = csv.DictWriter(f, result['colnames'])
            writer.writerow(result['values'])
            f.close()

        # Return
        return

    def models(self, models):
        """
        Set model.
        """
        # Copy models
        self._obs.models(models.copy())

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

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cstsdist(sys.argv)

    # Execute application
    app.execute()
