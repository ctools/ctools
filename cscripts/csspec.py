#! /usr/bin/env python
# ==========================================================================
# Spectral points generation script.
#
# Copyright (C) 2014-2015 Michael Mayer
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
    This class implements the creation of spectral points. It derives from
    the ctools.cscript class which provides support for parameter files,
    command line arguments, and logging. In that way the Python script
    behaves just as a regular ctool. 
    """
    def __init__(self, *argv):
        """
        Constructor.
        """

        # Set name
        self.name    = "csspec"
        self.version = "1.0.0"

        # Initialise some members
        self.obs  = gammalib.GObservations()
        self.fits = None

        # Initialise some members
        if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
            self.obs = argv[0]
            argv     = argv[1:]
        else:      
            self.obs = gammalib.GObservations()
            self.obs.clear()   
        self.m_outfile = ""

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
            pars.append(gammalib.GApplicationPar("inobs","f","a","events.fits","","","Input event list, counts cube, or observation definition XML file"))
            pars.append(gammalib.GApplicationPar("inmodel","f","a","$CTOOLS/share/models/crab.xml","","","Input model XML file"))
            pars.append(gammalib.GApplicationPar("srcname","s","a","Crab","","","Source name"))
            pars.append(gammalib.GApplicationPar("expcube","f","a","NONE","","","Input exposure cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("psfcube","f","a","NONE","","","Input PSF cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("bkgcube","f","a","NONE","","","Input background cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","prod2","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","South_0.5h","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("edisp","b","h","no","","","Apply energy dispersion?"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","spectrum.fits","","","Output spectrum file"))
            pars.append(gammalib.GApplicationPar("emin","r","a","0.1","","","Lower energy limit for spectral points (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","a","100.0","","","Upper energy limit for spectral points (TeV)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","20","1","10000","Number of spectral points"))
            pars.append(gammalib.GApplicationPar("ebinalg","s","a","LOG","FILE|LIN|LOG","", "Algorithm for defining energy bins"))
            pars.append(gammalib.GApplicationPar("ebinfile","f","a","NONE","","", "Name of the file containing the energy bin definition"))
            pars.append(gammalib.GApplicationPar("calc_ts","b","h","yes","yes|no","","Compute TS value in each bin"))
            pars.append(gammalib.GApplicationPar("calc_ulim","b","h","yes","yes|no","","Compute upper limit in each bin"))
            pars.append(gammalib.GApplicationPar("fix_srcs","b","h","yes","yes|no","","Fix other skymodel parameters"))
            pars.append(gammalib.GApplicationPar("fix_bkg","b","h","no","yes|no","","Fix background parameters"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csspec.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Set observation if not done before
        if self.obs == None or self.obs.size() == 0:
            self.require_inobs("csspec::get_parameters()")
            self.obs = self.get_observations()
        
        # Check if we have one binned cta observation, i.e. if we are in binned mode
        self.m_binned_mode = False    
        if self.obs.size() == 1:
            if self.obs[0].classname() == "GCTAObservation":
                if self.obs[0].eventtype() == "CountsCube":                
                    self.m_binned_mode = True        

        # Set models if we have none
        if self.obs.models().size() == 0:
            self.obs.models(self["inmodel"].filename())

        # Get source name   
        self.m_srcname = self["srcname"].string()

        # Set ebounds
        self.set_ebounds()
                 
        # Get edisp flag
        self.m_edisp = self["edisp"].boolean()

        # Read other parameters
        self.m_outfile = self["outfile"].filename()

        # Get other parameeters
        self.m_calc_ulimit = self["calc_ulim"].boolean()
        self.m_calc_ts     = self["calc_ts"].boolean()
        self.m_fix_bkg     = self["fix_bkg"].boolean()
        self.m_fix_srcs    = self["fix_srcs"].boolean()

        # Set some fixed parameters
        self.m_log     = False # Logging in client tools
        self.m_chatter = self["chatter"].integer()
        self.m_clobber = self["clobber"].boolean()
        self.m_debug   = self["debug"].boolean()

        # Return
        return

    def set_ebounds(self):
        """
        Set energy boundaries.
        """
        # Get ebounds
        if self.m_binned_mode:
            
            #  Extract cube ebounds
            cube_ebounds = self.obs[0].events().ebounds()
            
            # Read user parameters
            self.m_enumbins = self["enumbins"].integer()
            self.m_emin     = self["emin"].real() - 1e-6 # Rounding tolerance
            self.m_emax     = self["emax"].real() + 1e-6 # Rounding tolerance

            # Set energy range
            emin = gammalib.GEnergy(self.m_emin, "TeV")
            emax = gammalib.GEnergy(self.m_emax, "TeV")

            # Determine all cube layers that overlap with the [emin,emax]
            # energy range.
            use_layers = []
            for i in range(cube_ebounds.size()):
                if cube_ebounds.emin(i) < emax and cube_ebounds.emax(i) > emin:
                    use_layers.append(i)

            # If no layers overlap then throw an exception
            if (len(use_layers) == 0):
                msg = "Energy range ["+str(cube_ebounds.emin(i))+ \
                      ", "+str(cube_ebounds.emax(i))+"] of counts "+ \
                      "cube does not overlap with specified energy range ["+ \
                      str(emin)+", "+str(emax)+"]. Specify "+ \
                      "overlapping energy range."
                raise RuntimeError(msg)

            # Determine number of layers to use for each spectral bin
            n_layers = int(len(use_layers) / self.m_enumbins)

            # Don't allow zero layers
            if n_layers == 0:
                n_layers = 1

            # Create energy boundaries
            self.m_ebounds = gammalib.GEbounds()
            i_start = 0
            i_stop  = n_layers - 1
            while i_start < len(use_layers):
                if i_stop >= len(use_layers):
                    i_stop = len(use_layers) - 1
                ebinmin = cube_ebounds.emin(use_layers[i_start])
                ebinmax = cube_ebounds.emax(use_layers[i_stop])
                self.m_ebounds.append(ebinmin, ebinmax)
                i_start += n_layers
                i_stop  += n_layers
            #for i in range(len(use_layers))[::n_layers]:
            #    i_start = i
            #    i_stop  = i + n_layers - 1
            #    if i_stop >= len(use_layers):
            #        i_stop = len(use_layers) - 1
            #    ebinmin = cube_ebounds.emin(use_layers[i_start])
            #    ebinmax = cube_ebounds.emax(use_layers[i_stop])
            #    self.m_ebounds.append(ebinmin, ebinmax)
                    
        # Unbinned mode       
        else:
            self.m_ebounds = self.create_ebounds()

        # Return
        return

    def models(self, models):
        """
        Set model.
        """
        # Copy models
        self.obs.models(models.clone())

        # Return
        return

    def spectrum(self):
        """
        Return spectrum as FITS object.
        """
        # Return
        return self.fits

    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()

        # Save residual map
        self.fits.saveto(self.m_outfile, self.m_clobber)

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

        # Write spectral binning into header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Spectral binning")
            if self.m_binned_mode:
                cube_ebounds = self.obs[0].events().ebounds()
                self.log.parformat("Counts cube energy range")
                self.log(str(cube_ebounds.emin()))
                self.log(" - ")
                self.log(str(cube_ebounds.emax()))
                self.log("\n")
            for i in range(self.m_ebounds.size()):
                self.log.parformat("Bin "+str(i+1))
                self.log(str(self.m_ebounds.emin(i)))
                self.log(" - ")
                self.log(str(self.m_ebounds.emax(i)))
                self.log("\n")

        # Write observation into logger
        if self.logTerse():
            self.log("\n")
            self.log.header1("Observation")
            self.log(str(self.obs))
            self.log("\n")

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Adjust model parameters")

        # Adjust model parameters dependent on input user parameters
        for model in self.obs.models():

            # Set TS flag for all models to false.
            # Source of interest will be set to true later
            model.tscalc(False)

            # Log model name
            if self.logExplicit():
                self.log.header3(model.name())

            # Deal with the source of interest    
            if model.name() == self.m_srcname:
                for par in model:
                    if par.is_free() and self.logExplicit():
                        self.log(" Fixing \""+par.name()+"\"\n")
                    par.fix()
                normpar = model.spectral()[0]
                if normpar.is_fixed() and self.logExplicit():
                    self.log(" Freeing \""+normpar.name()+"\"\n")
                normpar.free()
                if self.m_calc_ts:
                    model.tscalc(True)

            elif self.m_fix_bkg and not model.classname() == "GModelSky":
                for par in model:
                    if par.is_free() and self.logExplicit():
                        self.log(" Fixing \""+par.name()+"\"\n")
                    par.fix()

            elif self.m_fix_srcs and model.classname() == "GModelSky":
                for par in model:
                    if par.is_free() and self.logExplicit():
                        self.log(" Fixing \""+par.name()+"\"\n")
                    par.fix()

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Generate spectrum")  
            self.log(str(self.m_ebounds))    

        # Initialise FITS Table with extension "SPECTRUM"
        table = gammalib.GFitsBinTable(self.m_ebounds.size())
        table.extname("SPECTRUM")

        # Add Header for compatibility with gammalib.GMWLSpectrum
        table.card("INSTRUME", "CTA", "Name of Instrument")
        table.card("TELESCOP", "CTA", "Name of Telescope")

        # Create FITS table columns
        energy       = gammalib.GFitsTableDoubleCol("Energy", self.m_ebounds.size())
        energy_low   = gammalib.GFitsTableDoubleCol("ed_Energy", self.m_ebounds.size())
        energy_high  = gammalib.GFitsTableDoubleCol("eu_Energy", self.m_ebounds.size())
        flux         = gammalib.GFitsTableDoubleCol("Flux", self.m_ebounds.size())
        flux_err     = gammalib.GFitsTableDoubleCol("e_Flux", self.m_ebounds.size())
        TSvalues     = gammalib.GFitsTableDoubleCol("TS", self.m_ebounds.size())
        ulim_values  = gammalib.GFitsTableDoubleCol("UpperLimit", self.m_ebounds.size())
        Npred_values = gammalib.GFitsTableDoubleCol("Npred", self.m_ebounds.size())
        energy.unit("TeV")
        energy_low.unit("TeV")
        energy_high.unit("TeV")
        flux.unit("erg/cm2/s")
        flux_err.unit("erg/cm2/s")
        ulim_values.unit("erg/cm2/s")

        # Loop over energy bins
        for i in range(self.m_ebounds.size()):

            # Log information
            if self.logExplicit():
                self.log("\n")
                self.log.header2("Energy bin "+str(i+1))

            # Get energy boundaries
            emin      = self.m_ebounds.emin(i)
            emax      = self.m_ebounds.emax(i)
            elogmean  = self.m_ebounds.elogmean(i)
            elogmean2 = elogmean.MeV() * elogmean.MeV()    

            # Store energy as TeV
            energy[i] = elogmean.TeV()

            # Store energy errors
            energy_low[i]  = (elogmean - emin).TeV()
            energy_high[i] = (emax - elogmean).TeV()

            # use ctselect for unbinned analysis
            if not self.m_binned_mode:
                
                # Log information
                if self.logExplicit():
                    self.log.header3("Selecting events")
    
                # Select events
                select = ctools.ctselect(self.obs)
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

            # use ctcubemask for binned analysis
            else:

                # Header
                if self.logExplicit():
                    self.log.header3("Filtering cube")

                # Select layers
                cubemask            = ctools.ctcubemask(self.obs)
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
            if self.logExplicit():
                self.log.header3("Performing fit")

            # Likelihood
            like          = ctools.ctlike(obs)
            like["edisp"] = self.m_edisp
            like.run()

            # Skip bin if no event was present
            if like.obs().logL() == 0.0:

                # Log information
                if self.logExplicit():
                    self.log("No event in this bin. ")
                    self.log("Likelihood is zero. ")
                    self.log("Bin is skipped.")

                # Set all values to 0
                flux[i]         = 0.0
                flux_err[i]     = 0.0
                TSvalues[i]     = 0.0
                ulim_values[i]  = 0.0
                Npred_values[i] = 0.0
                continue

            # Get results
            fitted_models = like.obs().models()
            source        = fitted_models[self.m_srcname]

            # Calculate Upper Limit            
            ulimit_value = -1.0
            if self.m_calc_ulimit:

                # Logging information
                if self.logExplicit():
                    self.log.header3("Computing upper limit")

                # Create upper limit object  
                ulimit = ctools.ctulimit(like.obs())
                ulimit["srcname"] = self.m_srcname
                ulimit["eref"]    = elogmean.TeV()

                # Try to run upper limit and catch exceptions
                try:
                    ulimit.run()
                    ulimit_value = ulimit.diff_ulimit()
                except:
                    if self.logExplicit():
                        self.log("Upper limit calculation failed.")
                    ulimit_value = -1.0

            # Get TS value
            TS = -1.0
            if self.m_calc_ts:
                TS = source.ts() 

            # Compute Npred value (only works for unbinned analysis)
            Npred = 0.0
            if not self.m_binned_mode:
                for observation in like.obs():
                    Npred += observation.npred(source)  

            # Get differential flux    
            fitted_flux = source.spectral().eval(elogmean,gammalib.GTime())

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
            if self.logTerse(): 
                self.log("\n")
                self.log.parformat("Bin "+str(i+1))
                self.log(str(flux[i]))
                self.log(" +/- ")
                self.log(str(flux_err[i]))
                if self.m_calc_ulimit and ulim_values[i] > 0.0:
                    self.log(" [< "+str(ulim_values[i])+"]")
                self.log(" erg/cm2/s")
                if self.m_calc_ts and TSvalues[i] > 0.0:
                    self.log(" (TS = "+str(TS)+")")

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
        self.fits = gammalib.GFits()
        self.fits.append(table)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates spectral points.
    """
    # Create instance of application
    app = csspec(sys.argv)

    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
