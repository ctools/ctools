#! /usr/bin/env python
# ==========================================================================
# Generates a spectrum.
#
# Copyright (C) 2014-2016 Michael Mayer
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


# ============ #
# csspec class #
# ============ #
class csspec(ctools.cscript):
    """
    Generates a spectrum.
    
    This class implements the creation of spectral points. It derives from
    the ctools.cscript class which provides support for parameter files,
    command line arguments, and logging. In that way the Python script
    behaves just as a regular ctool. 
    """

    # Constructors and destructors
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = "csspec"
        self._version = "1.1.0"

        # Initialise some members
        self._ebounds     = gammalib.GEbounds()
        self._fits        = None
        self._binned_mode = False
        self._srcname     = ""
        self._edisp       = False
        self._calc_ulimit = True
        self._calc_ts     = True
        self._fix_bkg     = False
        self._fix_srcs    = True
        self._chatter     = 2
        self._clobber     = True
        self._debug       = False

        # Initialise observation container from constructor arguments.
        self._obs = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return

    def __del__(self):
        """
        Destructor.
        """
        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Set observation if not done before
        if self._obs == None or self._obs.size() == 0:
            self._require_inobs("csspec::get_parameters()")
            self._obs = self._get_observations()
        
        # Check if we have one binned cta observation, i.e. if we are in binned mode
        self._binned_mode = False    
        if self._obs.size() == 1:
            if self._obs[0].classname() == "GCTAObservation":
                if self._obs[0].eventtype() == "CountsCube":                
                    self._binned_mode = True        

        # Set models if we have none
        if self._obs.models().size() == 0:
            self._obs.models(self["inmodel"].filename())

        # Get source name   
        self._srcname = self["srcname"].string()

        # Set ebounds
        self._set_ebounds()
                 
        # Get edisp flag
        self._edisp = self["edisp"].boolean()

        # Get other parameeters
        self._calc_ulimit = self["calc_ulim"].boolean()
        self._calc_ts     = self["calc_ts"].boolean()
        self._fix_bkg     = self["fix_bkg"].boolean()
        self._fix_srcs    = self["fix_srcs"].boolean()

        # Set some fixed parameters
        self._chatter = self["chatter"].integer()
        self._clobber = self["clobber"].boolean()
        self._debug   = self["debug"].boolean()

        # Read ahead output parameters
        if self._read_ahead():
            self["outfile"].filename()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Return
        return

    def _set_ebounds(self):
        """
        Set energy boundaries.
        """
        # Get ebounds
        if self._binned_mode:
            
            #  Extract cube ebounds
            cube_ebounds = self._obs[0].events().ebounds()
            
            # Read user parameters
            enumbins = self["enumbins"].integer()
            emin_val = self["emin"].real() - 1e-6 # Rounding tolerance
            emax_val = self["emax"].real() + 1e-6 # Rounding tolerance

            # Set energy range
            emin = gammalib.GEnergy(emin_val, "TeV")
            emax = gammalib.GEnergy(emax_val, "TeV")

            # Determine all cube layers that overlap with the [emin,emax]
            # energy range.
            use_layers = []
            for i in range(cube_ebounds.size()):
                if cube_ebounds.emin(i) < emax and cube_ebounds.emax(i) > emin:
                    use_layers.append(i)

            # If no layers overlap then throw an exception
            if (len(use_layers) == 0):
                msg = "Energy range ["+str(cube_ebounds.emin())+ \
                      ", "+str(cube_ebounds.emax())+"] of counts "+ \
                      "cube does not overlap with specified energy "+ \
                      "range ["+ \
                      str(emin)+", "+str(emax)+"]. Specify "+ \
                      "overlapping energy range."
                raise RuntimeError(msg)

            # Determine number of layers to use for each spectral bin
            n_layers = int(len(use_layers) / enumbins)

            # Don't allow zero layers
            if n_layers == 0:
                n_layers = 1

            # Create energy boundaries
            self._ebounds = gammalib.GEbounds()
            i_start = 0
            i_stop  = n_layers - 1
            while i_start < len(use_layers):
                if i_stop >= len(use_layers):
                    i_stop = len(use_layers) - 1
                ebinmin = cube_ebounds.emin(use_layers[i_start])
                ebinmax = cube_ebounds.emax(use_layers[i_stop])
                self._ebounds.append(ebinmin, ebinmax)
                i_start += n_layers
                i_stop  += n_layers
                    
        # Unbinned mode       
        else:
            self._ebounds = self._create_ebounds()

        # Return
        return


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

        # Write spectral binning into header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Spectral binning")
            if self._binned_mode:
                cube_ebounds = self._obs[0].events().ebounds()
                self._log.parformat("Counts cube energy range")
                self._log(str(cube_ebounds.emin()))
                self._log(" - ")
                self._log(str(cube_ebounds.emax()))
                self._log("\n")
            for i in range(self._ebounds.size()):
                self._log.parformat("Bin "+str(i+1))
                self._log(str(self._ebounds.emin(i)))
                self._log(" - ")
                self._log(str(self._ebounds.emax(i)))
                self._log("\n")

        # Write observation into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Observation")
            self._log(str(self._obs))
            self._log("\n")

        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Adjust model parameters")

        # Adjust model parameters dependent on input user parameters
        for model in self._obs.models():

            # Set TS flag for all models to false.
            # Source of interest will be set to true later
            model.tscalc(False)

            # Log model name
            if self._logExplicit():
                self._log.header3(model.name())

            # Deal with the source of interest    
            if model.name() == self._srcname:
                for par in model:
                    if par.is_free() and self._logExplicit():
                        self._log(" Fixing \""+par.name()+"\"\n")
                    par.fix()
                normpar = model.spectral()[0]
                if normpar.is_fixed() and self._logExplicit():
                    self._log(" Freeing \""+normpar.name()+"\"\n")
                normpar.free()
                if self._calc_ts:
                    model.tscalc(True)

            elif self._fix_bkg and not model.classname() == "GModelSky":
                for par in model:
                    if par.is_free() and self._logExplicit():
                        self._log(" Fixing \""+par.name()+"\"\n")
                    par.fix()

            elif self._fix_srcs and model.classname() == "GModelSky":
                for par in model:
                    if par.is_free() and self._logExplicit():
                        self._log(" Fixing \""+par.name()+"\"\n")
                    par.fix()

        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Generate spectrum")  
            self._log(str(self._ebounds))    

        # Initialise FITS Table with extension "SPECTRUM"
        table = gammalib.GFitsBinTable(self._ebounds.size())
        table.extname("SPECTRUM")

        # Add Header for compatibility with gammalib.GMWLSpectrum
        table.card("INSTRUME", "CTA", "Name of Instrument")
        table.card("TELESCOP", "CTA", "Name of Telescope")

        # Create FITS table columns
        nrows        = self._ebounds.size()
        energy       = gammalib.GFitsTableDoubleCol("Energy", nrows)
        energy_low   = gammalib.GFitsTableDoubleCol("ed_Energy", nrows)
        energy_high  = gammalib.GFitsTableDoubleCol("eu_Energy", nrows)
        flux         = gammalib.GFitsTableDoubleCol("Flux", nrows)
        flux_err     = gammalib.GFitsTableDoubleCol("e_Flux", nrows)
        TSvalues     = gammalib.GFitsTableDoubleCol("TS", nrows)
        ulim_values  = gammalib.GFitsTableDoubleCol("UpperLimit", nrows)
        Npred_values = gammalib.GFitsTableDoubleCol("Npred", nrows)
        energy.unit("TeV")
        energy_low.unit("TeV")
        energy_high.unit("TeV")
        flux.unit("erg/cm2/s")
        flux_err.unit("erg/cm2/s")
        ulim_values.unit("erg/cm2/s")

        # Loop over energy bins
        for i in range(nrows):

            # Log information
            if self._logExplicit():
                self._log("\n")
                self._log.header2("Energy bin "+str(i+1))

            # Get energy boundaries
            emin      = self._ebounds.emin(i)
            emax      = self._ebounds.emax(i)
            elogmean  = self._ebounds.elogmean(i)
            elogmean2 = elogmean.MeV() * elogmean.MeV()    

            # Store energy as TeV
            energy[i] = elogmean.TeV()

            # Store energy errors
            energy_low[i]  = (elogmean - emin).TeV()
            energy_high[i] = (emax - elogmean).TeV()

            # Use ctselect for unbinned analysis
            if not self._binned_mode:
                
                # Log information
                if self._logExplicit():
                    self._log.header3("Selecting events")
    
                # Select events
                select = ctools.ctselect(self._obs)
                select["emin"] = emin.TeV()    
                select["emax"] = emax.TeV() 
                select["tmin"] = "UNDEFINED"
                select["tmax"] = "UNDEFINED"
                select["rad"]  = "UNDEFINED"
                select["ra"]   = "UNDEFINED"
                select["dec"]  = "UNDEFINED"
                select.run()  
    
                # Retrieve observation
                obs = select.obs()

            # Use ctcubemask for binned analysis
            else:

                # Header
                if self._logExplicit():
                    self._log.header3("Filtering cube")

                # Select layers
                cubemask            = ctools.ctcubemask(self._obs)
                cubemask["regfile"] = "NONE"
                cubemask["ra"]      = "UNDEFINED"
                cubemask["dec"]     = "UNDEFINED"
                cubemask["rad"]     = "UNDEFINED"
                cubemask["emin"]    = emin.TeV() 
                cubemask["emax"]    = emax.TeV()
                cubemask.run() 
                
                # Set new binned observation
                obs = cubemask.obs()

            # Header
            if self._logExplicit():
                self._log.header3("Performing fit")

            # Likelihood
            like          = ctools.ctlike(obs)
            like["edisp"] = self._edisp
            like.run()

            # Skip bin if no event was present
            if like.obs().logL() == 0.0:

                # Log information
                if self._logExplicit():
                    self._log("No event in this bin. ")
                    self._log("Likelihood is zero. ")
                    self._log("Bin is skipped.")

                # Set all values to 0
                flux[i]         = 0.0
                flux_err[i]     = 0.0
                TSvalues[i]     = 0.0
                ulim_values[i]  = 0.0
                Npred_values[i] = 0.0
                continue

            # Get results
            fitted_models = like.obs().models()
            source        = fitted_models[self._srcname]

            # Calculate Upper Limit            
            ulimit_value = -1.0
            if self._calc_ulimit:

                # Logging information
                if self._logExplicit():
                    self._log.header3("Computing upper limit")

                # Create upper limit object  
                ulimit = ctools.ctulimit(like.obs())
                ulimit["srcname"] = self._srcname
                ulimit["eref"]    = elogmean.TeV()

                # Try to run upper limit and catch exceptions
                try:
                    ulimit.run()
                    ulimit_value = ulimit.diff_ulimit()
                except:
                    if self._logExplicit():
                        self._log("Upper limit calculation failed.")
                    ulimit_value = -1.0

            # Get TS value
            TS = -1.0
            if self._calc_ts:
                TS = source.ts() 

            # Compute Npred value (only works for unbinned analysis)
            Npred = 0.0
            if not self._binned_mode:
                for observation in like.obs():
                    Npred += observation.npred(source)  

            # Get differential flux    
            fitted_flux = source.spectral().eval(elogmean)

            # Compute flux error
            parvalue  = source.spectral()[0].value()
            rel_error = source.spectral()[0].error() / parvalue        
            e_flux    = fitted_flux * rel_error

            # Set values for storage
            TSvalues[i] = TS

            # Set npred values 
            Npred_values[i] = Npred

            # Convert fluxes to nuFnu
            flux[i]     = fitted_flux * elogmean2 * gammalib.MeV2erg
            flux_err[i] = e_flux      * elogmean2 * gammalib.MeV2erg
            if ulimit_value > 0.0:
                ulim_values[i] = ulimit_value * elogmean2 * gammalib.MeV2erg

            # Log information
            if self._logTerse(): 
                self._log("\n")
                self._log.parformat("Bin "+str(i+1))
                self._log(str(flux[i]))
                self._log(" +/- ")
                self._log(str(flux_err[i]))
                if self._calc_ulimit and ulim_values[i] > 0.0:
                    self._log(" [< "+str(ulim_values[i])+"]")
                self._log(" erg/cm2/s")
                if self._calc_ts and TSvalues[i] > 0.0:
                    self._log(" (TS = "+str(TS)+")")

        # Append filled columns to fits table    
        table.append(energy)
        table.append(energy_low)
        table.append(energy_high)
        table.append(flux)
        table.append(flux_err)
        table.append(TSvalues)
        table.append(ulim_values)
        table.append(Npred_values)

        # Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append(table)

        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Open logfile
        self._logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save spectrum
        self.save()

        # Return
        return

    def save(self):
        """
        Save spectrum.
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save spectrum")

        # Get outmap parameter
        outfile = self["outfile"].filename()
        
        # Continue only filename and residual map are valid
        if self._fits != None:

            # Log file name
            if self._logTerse():
                self._log(gammalib.parformat("Spectrum file"))
                self._log(outfile.url())
                self._log("\n")

            # Save spectrum
            self._fits.saveto(outfile, self._clobber)

        # Return
        return

    def spectrum(self):
        """
        Return spectrum FITS file.

        Returns:
            FITS file containing spectrum.
        """
        # Return
        return self._fits

    def models(self, models):
        """
        Set model.
        """
        # Copy models
        self._obs.models(models.clone())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csspec(sys.argv)

    # Execute application
    app.execute()
