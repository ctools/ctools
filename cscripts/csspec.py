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
        self.obs = None 

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
            pars.append(gammalib.GApplicationPar("inobs","f","a","events.fits","","","Event list, counts cube, or observation definition file"))
            pars.append(gammalib.GApplicationPar("inmodel","f","a","$CTOOLS/share/models/crab.xml","","","Source model"))
            pars.append(gammalib.GApplicationPar("srcname","s","a","Crab","","","Source name"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","spectrum.fits","","","Output file name"))
            pars.append(gammalib.GApplicationPar("expcube","f","a","NONE","","","Exposure cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("psfcube","f","a","NONE","","","PSF cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("bkgcube","s","a","NONE","","","Background cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","prod2","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","South_50h","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("emin","r","h","0.1","","","Lower energy limit for spectral points(TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","h","100.0","","","Upper energy limit for spectral points(TeV)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","20","","","Number of spectral points"))
            pars.append(gammalib.GApplicationPar("ebinalg","s", "h","LOG","FILE|LIN|LOG","", "Algorithm for defining energy bins"))
            pars.append(gammalib.GApplicationPar("binned","b","a","no","yes|no","","Use binned analysis in each energy bin"))
            pars.append(gammalib.GApplicationPar("nebins","i","h","5","","","Number of energy bins per spectral point"))
            pars.append(gammalib.GApplicationPar("coordsys","s","a","CEL","CEL|GAL","","Coordinate System"))
            pars.append(gammalib.GApplicationPar("proj","s","a","CAR","AIT|AZP|CAR|MER|MOL|STG|TAN","","Projection method"))
            pars.append(gammalib.GApplicationPar("xref","r","a","83.63","0","360","First coordinate of image center in degrees (RA or galactic l)"))
            pars.append(gammalib.GApplicationPar("yref","r","a","22.01","-90","90","Second coordinate of image center in degrees (DEC or galactic b)"))
            pars.append(gammalib.GApplicationPar("nxpix","i","a","200","","","Size of the X axis in pixels"))
            pars.append(gammalib.GApplicationPar("nypix","i","a","200","","","Size of the Y axis in pixels"))
            pars.append(gammalib.GApplicationPar("binsz","r","a","0.02","","","Pixel size (deg/pixel)"))
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

        # Set models if we have none
        if self.obs.models().size() == 0:
            self.obs.models(self["inmodel"].filename())

        # Get source name   
        self.m_srcname = self["srcname"].string()

        # Get ebounds             
        self.m_ebounds = self.create_ebounds()

        # Get binning flag
        self.m_binned = self["binned"].boolean()
        if self.m_binned:
            self.m_xref     = self["xref"].real()
            self.m_yref     = self["yref"].real()
            self.m_nxpix    = self["nxpix"].integer()
            self.m_nypix    = self["nypix"].integer()
            self.m_binsz    = self["binsz"].real()
            self.m_coordsys = self["coordsys"].string()
            self.m_proj     = self["proj"].string()
            self.m_ebins    = self["nebins"].integer()

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

    def models(self, models):
        """
        Set model.
        """
        # Copy models
        self.obs.models(models.clone())

        # Return
        return

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
                if par.is_fixed() and self.logExplicit():
                    self.log(" Freeing \""+par.name()+"\"\n")
                model.spectral()[0].free()
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

        # Initialise FITS Table with extension "SPECTRUM"
        table = gammalib.GFitsBinTable(self.m_ebounds.size())
        table.extname("SPECTRUM")

        # Add Header for compatibility with gammalib.GMWLSpectrum
        table.card("INSTRUME", "CTA", "Name of Instrument")
        table.card("TELESCOP", "CTA", "Name of Telescope")

        # Create FITS table columns
        energy       = gammalib.GFitsTableDoubleCol("Energy", self.m_ebounds.size())
        energy.unit("TeV")
        energy_low   = gammalib.GFitsTableDoubleCol("ed_Energy", self.m_ebounds.size())
        energy_low.unit("TeV")
        energy_high  = gammalib.GFitsTableDoubleCol("eu_Energy", self.m_ebounds.size())
        energy_high.unit("TeV")
        flux         = gammalib.GFitsTableDoubleCol("Flux", self.m_ebounds.size())
        flux.unit("erg/cm2/s")
        flux_err     = gammalib.GFitsTableDoubleCol("e_Flux", self.m_ebounds.size())
        flux_err.unit("erg/cm2/s")
        TSvalues     = gammalib.GFitsTableDoubleCol("TS", self.m_ebounds.size())
        ulim_values  = gammalib.GFitsTableDoubleCol("UpperLimit", self.m_ebounds.size())
        ulim_values.unit("erg/cm2/s")
        Npred_values = gammalib.GFitsTableDoubleCol("Npred", self.m_ebounds.size())


        # Loop over energy bins
        for i in range(self.m_ebounds.size()):

            # Log information
            if self.logTerse():
                self.log("\n")
                self.log.header2("Energy bin "+str(i))

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

            # Binned analysis
            if self.m_binned:

                # Header
                if self.logTerse():
                    self.log.header3("Binning events")

                # Bin events
                bin = ctools.ctbin(select.obs())
                bin["usepnt"]   = False
                bin["ebinalg"]  = "LOG"
                bin["xref"]     = self.m_xref
                bin["yref"]     = self.m_yref
                bin["binsz"]    = self.m_binsz
                bin["nxpix"]    = self.m_nxpix
                bin["nypix"]    = self.m_nypix
                bin["enumbins"] = self.m_ebins
                bin["emin"]     = emin.TeV()
                bin["emax"]     = emax.TeV()        
                bin["coordsys"] = self.m_coordsys
                bin["proj"]     = self.m_proj
                bin.run()

                # Header
                if self.logTerse():
                    self.log.header3("Creating exposure cube")

                # Create exposure cube
                expcube = ctools.ctexpcube(select.obs())
                expcube["incube"]   = "NONE"
                expcube["usepnt"]   = False
                expcube["ebinalg"]  = "LOG"
                expcube["xref"]     = self.m_xref
                expcube["yref"]     = self.m_yref
                expcube["binsz"]    = self.m_binsz
                expcube["nxpix"]    = self.m_nxpix
                expcube["nypix"]    = self.m_nypix
                expcube["enumbins"] = self.m_ebins
                expcube["emin"]     = emin.TeV()
                expcube["emax"]     = emax.TeV() 
                expcube["coordsys"] = self.m_coordsys
                expcube["proj"]     = self.m_proj               
                expcube.run()

                # Header
                if self.logTerse():
                    self.log.header3("Creating PSF cube")

                # Create psf cube
                psfcube = ctools.ctpsfcube(select.obs())
                psfcube["incube"]   = "NONE"
                psfcube["usepnt"]   = False
                psfcube["ebinalg"]  = "LOG"
                psfcube["xref"]     = self.m_xref
                psfcube["yref"]     = self.m_yref
                psfcube["binsz"]    = self.m_binsz
                psfcube["nxpix"]    = self.m_nxpix
                psfcube["nypix"]    = self.m_nypix
                psfcube["enumbins"] = self.m_ebins
                psfcube["emin"]     = emin.TeV()
                psfcube["emax"]     = emax.TeV()  
                psfcube["coordsys"] = self.m_coordsys
                psfcube["proj"]     = self.m_proj               
                psfcube.run()

                # Header
                if self.logTerse():
                    self.log.header3("Creating background cube")

                # Create background cube
                bkgcube = ctools.ctbkgcube(select.obs())
                bkgcube["incube"]   = "NONE"
                bkgcube["usepnt"]   = False
                bkgcube["ebinalg"]  = "LOG"
                bkgcube["xref"]     = self.m_xref
                bkgcube["yref"]     = self.m_yref
                bkgcube["binsz"]    = self.m_binsz
                bkgcube["nxpix"]    = self.m_nxpix
                bkgcube["nypix"]    = self.m_nypix
                bkgcube["enumbins"] = self.m_ebins
                bkgcube["emin"]     = emin.TeV()
                bkgcube["emax"]     = emax.TeV() 
                bkgcube["coordsys"] = self.m_coordsys
                bkgcube["proj"]     = self.m_proj                
                bkgcube.run()

                # Set new binned observation
                obs = bin.obs()

                # Set precomputed binned response
                obs[0].response(expcube.expcube(), psfcube.psfcube(), bkgcube.bkgcube())

                # Get new models
                models = bkgcube.models()

                # Fix background models if required
                if self.m_fix_bkg:
                    for model in models:
                        if not model.classname() == "GModelSky":
                            for par in model:
                                par.fix()

                # Set new models to binned observation           
                obs.models(models)

            # Header
            if self.logTerse():
                self.log.header3("Performing fit")

            # Likelihood
            like = ctools.ctlike(obs)
            like.run()

            # Skip bin if no event was present
            if like.obs().logL() == 0.0:

                # Log information
                if self.logTerse():
                    self.log("No event in this bin. Bin is skipped\n")

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
                if self.logTerse():
                    self.log.header3("Computing upper limit")

                # Create upper limit object  
                ulimit = ctools.ctulimit(like.obs())
                ulimit["srcname"] = self.m_srcname
                ulimit["eref"] = elogmean.TeV()

                # Try to run upper limit and catch exceptions
                try:
                    ulimit.run()
                    ulimit_value = ulimit.diff_ulimit()
                except:
                    if self.logTerse():
                        self.log("Upper limit calculation failed\n")
                    ulimit_value = -1.0

            # Get TS value
            TS = -1.0
            if self.m_calc_ts:
                TS = source.ts() 

            # Compute Npred value
            Npred = 0.0
            for observation in like.obs():
                Npred += observation.npred(source)  

            # Get differential flux    
            fitted_flux = source.spectral().eval(elogmean,gammalib.GTime())

            # Compute flux error
            parvalue  = source.spectral()[0].value()
            rel_error = source.spectral()[0].error()/parvalue        
            e_flux    = fitted_flux*rel_error

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
            if self.logExplicit(): 
                self.log("Bin "+str(i)+" ["+str(emin.TeV())+"-"+str(emax.TeV())+"] TeV: ")
                self.log("Flux = "+str(flux[i]))
                self.log(" +- "+str(flux_err[i])+" [erg/cm2/s]")
                if self.m_calc_ts and TSvalues[i] > 0.0:
                    self.log(", TS = "+str(TS))
                if self.m_calc_ulimit and ulim_values[i] > 0.0:
                    self.log(", UL = "+str(ulim_values[i])+" [erg/cm2/s]")
                self.log("\n")

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
