#! /usr/bin/env python
# ==========================================================================
# This script generates the pull distribution for all model parameters.
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


# ============ #
# cspull class #
# ============ #
class cspull(ctools.cscript):
    """
    Generates pull distributions for a model.
    
    This class implements the pull distribution generation script. It derives
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
        self._name    = "cspull"
        self._version = "1.1.0"

        # Initialise some members
        self._edisp       = False
        self._exposure    = None
        self._psfcube     = None
        self._bckcube     = None
        self._edispcube   = None
        self._stackmodels = None
        self._coordsys    = "CEL"
        self._proj        = "TAN"
        self._log_clients = False
        self._offset      = 0.0
        self._ntrials     = 0
        self._outfile     = gammalib.GFilename()
        self._npix        = 0
        self._binsz       = 0.0
        self._pattern     = "single"
        self._enumbins    = 0
        self._chatter     = 2

        # Initialise observation container from constructor arguments.
        self._obs, argv = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile.
        """
        # If there are no observations in container then get some ...
        if self._obs.size() == 0:
            self._obs = self._get_observations()

            # Check for requested pattern and use above
            # observation parameters to set wobble pattern
            self._pattern = self["pattern"].string()
            if self._pattern == "four":
                self._obs = self._set_obs()

        # ... otherwise add response information and energy boundaries
        # in case they are missing.
        else:
            self._setup_observations(self._obs)

        # Get number of energy bins
        self._enumbins = self["enumbins"].integer()

        # Read parameters for binned if requested
        if self._enumbins != 0:
            self._npix     = self["npix"].integer()
            self._binsz    = self["binsz"].real()
            self._coordsys = self["coordsys"].string()
            self._proj     = self["proj"].string()

        # Set models if we have none
        if self._obs.models().size() == 0:
            self._obs.models(self["inmodel"].filename())

        # Read other parameters    
        self._outfile = self["outfile"].filename()
        self._ntrials = self["ntrials"].integer()
        self._edisp   = self["edisp"].boolean()
        self._offset  = self["offset"].real()
        self._chatter = self["chatter"].integer()

        # Query some parameters
        self["profile"].boolean()

        #  Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Return
        return

    def _set_obs(self):
        """
        Set CTA observation(s).

        Returns:
            An observation container.
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

    def _set_stacked_irf(self):
        """
        For stacked analysis prepare stacked irfs
        """
        # Write header into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Compute stacked response")

        # Get stacked exposure
        expcube = ctools.ctexpcube(self._obs)
        expcube["incube"]   = "NONE"
        expcube["usepnt"]   = True
        expcube["ebinalg"]  = "LOG"
        expcube["binsz"]    = self._binsz
        expcube["nxpix"]    = self._npix
        expcube["nypix"]    = self._npix
        expcube["enumbins"] = self._enumbins
        expcube["emin"]     = self["emin"].real()
        expcube["emax"]     = self["emax"].real()
        expcube["coordsys"] = self._coordsys
        expcube["proj"]     = self._proj
        expcube.run()

        # Notify exposure computation
        if self._logTerse():
            self._log("Computed exposure cube\n")

        # Get stacked Psf
        psfcube = ctools.ctpsfcube(self._obs)
        psfcube["incube"]   = "NONE"
        psfcube["usepnt"]   = True
        psfcube["ebinalg"]  = "LOG"
        psfcube["binsz"]    = self._binsz*10.0
        psfcube["nxpix"]    = self._npix/10
        psfcube["nypix"]    = self._npix/10
        psfcube["enumbins"] = self._enumbins
        psfcube["emin"]     = self["emin"].real()
        psfcube["emax"]     = self["emax"].real()
        psfcube["coordsys"] = self._coordsys
        psfcube["proj"]     = self._proj
        psfcube.run()

        # Notify Psf computation
        if self._logTerse():
            self._log("Computed point spread function cube\n")

        # Optionally get stacked Edisp
        if self._edisp:
            edispcube = ctools.ctedispcube(self._obs)
            edispcube["incube"]   = "NONE"
            edispcube["usepnt"]   = True
            edispcube["ebinalg"]  = "LOG"
            edispcube["binsz"]    = self._binsz*10.0
            edispcube["nxpix"]    = self._npix/10
            edispcube["nypix"]    = self._npix/10
            edispcube["enumbins"] = self._enumbins
            edispcube["emin"]     = self["emin"].real()
            edispcube["emax"]     = self["emax"].real()
            edispcube["coordsys"] = self._coordsys
            edispcube["proj"]     = self._proj
            edispcube.run()   

            # Store result
            self._edispcube = edispcube.edispcube().copy()
            
            # Logging
            if self._logTerse():
                self._log("Computed energy dispersion cube\n")

        # Get stacked background
        bkgcube = ctools.ctbkgcube(self._obs)
        bkgcube["incube"]   = "NONE"
        bkgcube["usepnt"]   = True
        bkgcube["ebinalg"]  = "LOG"
        bkgcube["binsz"]    = self._binsz
        bkgcube["nxpix"]    = self._npix
        bkgcube["nypix"]    = self._npix
        bkgcube["enumbins"] = self._enumbins
        bkgcube["emin"]     = self["emin"].real()
        bkgcube["emax"]     = self["emax"].real()
        bkgcube["coordsys"] = self._coordsys
        bkgcube["proj"]     = self._proj
        bkgcube.run()

        # Notify background cube computation
        if self._logTerse():
            self._log("Computed background cube\n")

        # Store results
        self._exposure    = expcube.expcube().copy()
        self._psfcube     = psfcube.psfcube().copy()
        self._bckcube     = bkgcube.bkgcube().copy()
        self._stackmodels = bkgcube.models().copy()

        # Return
        return

    def _trial(self, seed):
        """
        Compute the pull for a single trial.

        Args:
            seed: Random number generator seed.
        """
        # Write header
        if self._logNormal():
            self._log.header2("Trial "+str(seed+1))

        # Simulate events
        obs = obsutils.sim(self._obs,
                           nbins=self._enumbins,
                           seed=seed,
                           binsz=self._binsz,
                           npix=self._npix,
                           proj=self._proj,
                           coord=self._coordsys,
                           edisp=self._edisp,
                           log=self._log_clients,
                           debug=self._logDebug(),
                           chatter=self._chatter)

        # Set response for a stacked observation
        if self._obs.size() > 1 and self._enumbins > 0:
            if self._edisp:
                obs[0].response(self._exposure,  self._psfcube,
                                self._edispcube, self._bckcube)
            else:
                obs[0].response(self._exposure, self._psfcube,
                                self._bckcube)
            obs.models(self._stackmodels)

        # Determine number of events in simulation
        nevents = 0.0
        for run in obs:
            nevents += run.events().number()

        # Write simulation results
        if self._logNormal():
            self._log.header3("Simulation")
            self._log.parformat("Number of simulated events")
            self._log(nevents)
            self._log("\n")

        # Fit model
        if self["profile"].boolean():
            models = self._obs.models()
            for i in range(models.size()):
                model_name = models[i].name()
                like       = obsutils.cterror(self._obs, model_name,
                                              log=self._log_clients,
                                              debug=self._logDebug(),
                                              chatter=self._chatter)
        else:
            like = obsutils.fit(obs, edisp=self._edisp,
                                log=self._log_clients,
                                debug=self._logDebug(),
                                chatter=self._chatter)

        # Store results
        logL   = like.opt().value()
        npred  = like.obs().npred()
        models = like.obs().models()

        # Write result header
        if self._logNormal():
            self._log.header3("Pulls")

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
                    real_value   = self._obs.models()[i][k].value()
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
                    if self._logNormal():
                        self._log.parformat(name)
                        self._log(pull)
                        self._log(" (")
                        self._log(fitted_value)
                        self._log(" +/- ")
                        self._log(error)
                        self._log(")\n")

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

        # Write observation into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Observation(s)")
            self._log(str(self._obs))
            self._log("\n")

        # If several observations and binned: prepare stacked irfs
        if self._obs.size() > 1 and self._enumbins > 0:
            self._set_stacked_irf()

        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Generate pull distribution")

        # Loop over trials
        for seed in range(self._ntrials):

            # Make a trial
            result = self._trial(seed)

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

    def execute(self):
        """
        Execute the script.
        """
        # Open logfile
        self.logFileOpen()

        # Run the script
        self.run()

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


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cspull(sys.argv)

    # Execute application
    app.execute()
