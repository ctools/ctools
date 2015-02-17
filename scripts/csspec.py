#! /usr/bin/env python
# ==========================================================================
# spectral points generation script.
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

# ============== #
# cstsdist class #
# ============== #
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
        self.obs       = None 
              
        # Initialise some members
        if isinstance(argv[0],gammalib.GObservations):
            self.obs = argv[0]
            argv = argv[1:]
        else:      
            self.obs      = gammalib.GObservations()
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
            pars.append(gammalib.GApplicationPar("caldb","s","a","dummy","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","cta_dummy_irf","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("emin","r","h","0.1","","","Lower energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","h","100.0","","","Upper energy limit (TeV)"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","20","","","Number of energy bins"))
            pars.append(gammalib.GApplicationPar("ebinalg","s","h","LOG","LIN|LOG|FILE","","Binning algorithm"))
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
        # Get parameters
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

        # Read other parameters
        self.m_outfile   = self["outfile"].filename()
        
        # Get other parameeters
        self.m_calc_ulimit = self["calc_ulim"].boolean()
        self.m_calc_ts = self["calc_ts"].boolean()
        self.m_fix_bkg = self["fix_bkg"].boolean()
        self.m_fix_srcs = self["fix_srcs"].boolean()
             
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

        for model in self.obs.models():
            model.tscalc(False)
            if model.name() == self.m_srcname:
                for par in model:
                    par.fix()
                model.spectral()[0].free()
                if self.m_calc_ts:
                    model.tscalc(True)
                
            elif self.m_fix_bkg and not model.classname() == "GModelData":
                for par in model:
                    par.fix()
        
            elif self.m_fix_srcs and model.classname() == "GModelSky":
                for par in model:
                    par.fix()
        
        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Generate spectrum")      
        
        # Initialise FITS Table with extension "SPECTRUM"
        table = gammalib.GFitsBinTable(self.m_ebounds.size())
        table.extname("SPECTRUM")
        
        # Add Header for compatibility with gammalib.GMWLSpectrum
        table.card("Instrument","CTA","Name of Instrument")
        table.card("Telescope","CTA","Name of Telescope")
        
        # Create FITS table columns
        energy = gammalib.GFitsTableDoubleCol("energy",self.m_ebounds.size())
        energy.unit("TeV")
        energy_err = gammalib.GFitsTableDoubleCol("energy_err",self.m_ebounds.size())
        energy_err.unit("TeV")
        flux = gammalib.GFitsTableDoubleCol("flux",self.m_ebounds.size())
        flux.unit("erg/cm2/s")
        flux_err = gammalib.GFitsTableDoubleCol("flux_err",self.m_ebounds.size())
        flux_err.unit("erg/cm2/s")
        TSvalues = gammalib.GFitsTableDoubleCol("TS",self.m_ebounds.size())
        ulim_values  = gammalib.GFitsTableDoubleCol("ulimit",self.m_ebounds.size())
        ulim_values.unit("erg/cm2/s")
        
        # MeV2erg
        MeV2erg = 1.60217e-6
            
        for i in range(self.m_ebounds.size()):
            
            emin = self.m_ebounds.emin(i)
            emax = self.m_ebounds.emax(i)
            emean = self.m_ebounds.emean(i)
            elogmean = self.m_ebounds.elogmean(i)
            elogmean2 = elogmean.MeV() * elogmean.MeV()
            
            # Select events
            select = ctools.ctselect(self.obs)
            select["emin"].real(emin.TeV())    
            select["emax"].real(emax.TeV()) 
            select["tmin"].value("UNDEFINED")
            select["tmax"].value("UNDEFINED")
            select["rad"].value("UNDEFINED")
            select["ra"].value("UNDEFINED")
            select["dec"].value("UNDEFINED")
            select.run()        
            
            # likelihood
            like = ctools.ctlike(select.obs())
            like.run()
            
            # Get results
            fitted_models = like.obs().models()
            source = fitted_models[self.m_srcname]
            
            # Calculate Upper Limit            
            ulimit_value = -1.0
            #if self.m_calc_ulimit:
            #    ulimit = ctools.ctulimit(like.obs())
            #    ulimit["eref"].real(elogmean.TeV())
            #    ulimit.run()
            #    ulimit_value.ulimit.diff_ulimit()
            
            # Get TS value
            TS = -1.0
            if self.m_calc_ts:
                TS = source.ts() 
              
            # Get differential flux    
            fitted_flux = source.spectral().eval(elogmean,gammalib.GTime())
            
            # Compute flux error
            parvalue = source.spectral()[0].value()
            rel_error = source.spectral()[0].error()/parvalue        
            e_flux = fitted_flux*rel_error
            
            # Set values for storage
            TSvalues[i] = TS
            
            # Store energy as TeV
            energy[i] = elogmean.TeV()
            
            # Store energy error
            energy_err[i] = (emax - emean).TeV()
            
            # Convert fluxes to nuFnu
            flux[i] = fitted_flux * elogmean2 * MeV2erg
            flux_err[i] = e_flux * elogmean2 * MeV2erg
            if ulimit_value > 0.0:
                ulim_values[i] = ulimit_value * elogmean2 * MeV2erg

            
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
            
        # Append filles columns to fits table    
        table.append(energy)
        table.append(energy_err)
        table.append(flux)
        table.append(flux_err)
        table.append(TSvalues)
        table.append(ulim_values)
        
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
    