#! /usr/bin/env python
# ==========================================================================
# Light curve generation script.
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
class cslightcrv(ctools.cscript):
    """
    This class implements the creation of a light curve. It derives from
    the ctools.cscript class which provides support for parameter files,
    command line arguments, and logging. In that way the Python script
    behaves just as a regular ctool. 
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        
        # Set name
        self.name    = "cslightcrv"
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
            pars.append(gammalib.GApplicationPar("outfile","f","a","lightcurve.fits","","","Output file name"))
            pars.append(gammalib.GApplicationPar("expcube","f","a","NONE","","","Exposure cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("psfcube","f","a","NONE","","","PSF cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("bkgcube","s","a","NONE","","","Background cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","prod2","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","South_50h","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("tbinalg","s","h","GTI","FILE|LIN|GTI","", "Algorithm for defining time bins"))
            pars.append(gammalib.GApplicationPar("tmin","r","a","51544.5","","days","Lightcurve start time [MJD]"))
            pars.append(gammalib.GApplicationPar("tmax","r","a","51544.6","","days","Lightcurve stop time [MJD]"))
            pars.append(gammalib.GApplicationPar("tbins","i","a","5","","","Number of time bins"))
            pars.append(gammalib.GApplicationPar("tbinfile","f","a","tbins.fits","","", "File defining the time binning"))
            #pars.append(gammalib.GApplicationPar("binned","b","a","no","yes|no","","Use binned analysis in each time bin"))
            pars.append(gammalib.GApplicationPar("enumbins","i","a","0","","","Number of energy bins per light curve bin (0=unbinned)"))
            pars.append(gammalib.GApplicationPar("emin","r","a","0.1","","","Lower energy limit of events (TeV)"))
            pars.append(gammalib.GApplicationPar("emax","r","a","100.0","","","Upper energy limit of events (TeV)"))
            pars.append(gammalib.GApplicationPar("coordsys","s","a","CEL","CEL|GAL","","Coordinate System"))
            pars.append(gammalib.GApplicationPar("proj","s","a","CAR","AIT|AZP|CAR|MER|MOL|STG|TAN","","Projection method"))
            pars.append(gammalib.GApplicationPar("xref","r","a","83.63","0","360","First coordinate of image center in degrees (RA or galactic l)"))
            pars.append(gammalib.GApplicationPar("yref","r","a","22.01","-90","90","Second coordinate of image center in degrees (DEC or galactic b)"))
            pars.append(gammalib.GApplicationPar("nxpix","i","a","200","","","Size of the X axis in pixels"))
            pars.append(gammalib.GApplicationPar("nypix","i","a","200","","","Size of the Y axis in pixels"))
            pars.append(gammalib.GApplicationPar("binsz","r","a","0.02","","","Pixel size (deg/pixel)"))
            pars.append(gammalib.GApplicationPar("calc_ts","b","h","yes","yes|no","","Compute TS value for each bin?"))
            pars.append(gammalib.GApplicationPar("calc_ulim","b","h","yes","yes|no","","Compute upper limit for each bin?"))
            pars.append(gammalib.GApplicationPar("fix_srcs","b","h","yes","yes|no","","Fix other sky model parameters?"))
            pars.append(gammalib.GApplicationPar("fix_bkg","b","h","no","yes|no","","Fix background model parameters?"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","cslightcrv.log","","","Log filename"))
            pars.save(parfile)
        
        # Return
        return

    def create_tbounds(self):
        
        # Create time bin container
        self.m_tbins = gammalib.GGti()

        # check for binning algorithm
        if self["tbinalg"].string() == "FILE":
            try:
                self.m_tbins.load(self["tbinfile"].filename())
            except:               
                csv = gammalib.GCsv(self["tbinfile"].filename())
                for i in range(csv.nrows()):
                    tmin = gammalib.GTime()
                    tmax = gammalib.GTime()
                    tmin.mjd(csv.real(i,0))
                    tmax.mjd(csv.real(i,1))
                    self.m_tbins.append(tmin,tmax)
        elif self["tbinalg"].string() == "LIN":
            # Use linear time binning
            time_min = self["tmin"].real()
            time_max = self["tmax"].real()
            nbins    = self["tbins"].integer()            
            step     = (time_max - time_min) / float(nbins)
            for i in range(nbins):
                tmin = gammalib.GTime()
                tmin.mjd(time_min + i*step)# ref)
                tmax = gammalib.GTime()
                tmax.mjd(time_min + (i+1)*step)#ref)
                self.m_tbins.append(tmin,tmax)
        elif self["tbinalg"].string() == "GTI":
            # Use GTIs of observations
            for obs in self.obs:
                for i in range(obs.events().gti().size()):
                    self.m_tbins.append(obs.events().gti().tstart(i),obs.events().gti().tstop(i))
        else:
            raise AttributeError("tbinalg=\""+self["tbinalg"].string()+"\" unkown. Must be one of \"FILE\", \"LIN\" or \"GTI\"")

        # Return
        return
  
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Set observation if not done before
        if self.obs == None or self.obs.size() == 0:
            self.require_inobs("cslightcrv::get_parameters()")
            self.obs = self.get_observations()
        
        # Set models if we have none
        if self.obs.models().size() == 0:
            self.obs.models(self["inmodel"].filename())
          
        # Get source name   
        self.m_srcname = self["srcname"].string()
        
        # Get time boundaries             
        self.create_tbounds()

        # Unbinned or binned analysis?
        self.m_ebins    = self["enumbins"].integer()
        if self.m_ebins == 0:
            self.m_binned = False
        else:
            self.m_binned = True

        # Get energy range
        self.m_emin = self["emin"].real()
        self.m_emax = self["emax"].real()

        # Get binning flag
        #self.m_binned = self["binned"].boolean()
        if self.m_binned:
            self.m_coordsys = self["coordsys"].string()
            self.m_proj     = self["proj"].string()
            self.m_xref     = self["xref"].real()
            self.m_yref     = self["yref"].real()
            self.m_nxpix    = self["nxpix"].integer()
            self.m_nypix    = self["nypix"].integer()
            self.m_binsz    = self["binsz"].real()
            #self.m_ebins    = self["enumbins"].integer()
            #self.m_emin     = self["emin"].real()
            #self.m_emax     = self["emax"].real()

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
            self.log.header1("Generate lightcurve")      
        
        # Initialise FITS Table with extension "LIGHTCURVE"
        table = gammalib.GFitsBinTable(self.m_tbins.size())
        table.extname("LIGHTCURVE")
        
        # Add Header for compatibility with gammalib.GMWLSpectrum
        table.card("INSTRUME", "CTA", "Name of Instrument")
        table.card("TELESCOP", "CTA", "Name of Telescope")
             
        # Create FITS table columns        
        MJD = gammalib.GFitsTableDoubleCol("MJD", self.m_tbins.size())
        MJD.unit("days")
        e_MJD = gammalib.GFitsTableDoubleCol("e_MJD", self.m_tbins.size())
        e_MJD.unit("days")
        
        # Create a FITS column for every free parameter
        columns = []
        for par in self.obs.models()[self.m_srcname]:
            if par.is_free():
                col = gammalib.GFitsTableDoubleCol(par.name(), self.m_tbins.size())
                col.unit(par.unit())
                columns.append(col)
                e_col = gammalib.GFitsTableDoubleCol("e_"+par.name(), self.m_tbins.size())
                e_col.unit(par.unit())
                columns.append(e_col)
        
        # Create TS and upper limit columns
        TSvalues    = gammalib.GFitsTableDoubleCol("TS", self.m_tbins.size())
        ulim_values = gammalib.GFitsTableDoubleCol("UpperLimit", self.m_tbins.size())
        ulim_values.unit("ph/cm2/s")

        # Loop over energy bins
        for i in range(self.m_tbins.size()):

            # Log information
            if self.logTerse():
                self.log("\n")
                self.log.header2("Time bin "+str(i))

            # Get time boundaries
            tmin = self.m_tbins.tstart(i)
            tmax = self.m_tbins.tstop(i)
            
            # Compute time bin center and time width
            tmean   = (tmin + tmax)
            tmean  *= 0.5
            twidth  = (tmax - tmin)
            twidth *= 0.5 

            # Store time as MJD
            MJD[i] = tmean.mjd()
            e_MJD[i] = twidth.days()
            
            # Log information
            if self.logExplicit():
                self.log.header3("Selecting events")
                     
            # Select events
            select = ctools.ctselect(self.obs)
            select["emin"].real(self.m_emin)    
            select["emax"].real(self.m_emax) 
            select["tmin"].real(tmin.convert(select.time_reference()))
            select["tmax"].real(tmax.convert(select.time_reference()))
            select["rad"].value("UNDEFINED")
            select["ra"].value("UNDEFINED")
            select["dec"].value("UNDEFINED")
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
                bin["usepnt"].boolean(False)
                bin["ebinalg"].string("LOG")
                bin["xref"].real(self.m_xref)
                bin["yref"].real(self.m_yref)
                bin["binsz"].real(self.m_binsz)
                bin["nxpix"].integer(self.m_nxpix)
                bin["nypix"].integer(self.m_nypix)
                bin["enumbins"].integer(self.m_ebins)
                bin["emin"].real(self.m_emin)
                bin["emax"].real(self.m_emax)        
                bin["coordsys"].string(self.m_coordsys)
                bin["proj"].string(self.m_proj)            
                bin.run()
                
                # Header
                if self.logTerse():
                    self.log.header3("Creating exposure cube")
                
                # Create exposure cube
                expcube = ctools.ctexpcube(select.obs())
                expcube["incube"].filename("NONE")
                expcube["usepnt"].boolean(False)
                expcube["ebinalg"].string("LOG")
                expcube["xref"].real(self.m_xref)
                expcube["yref"].real(self.m_yref)
                expcube["binsz"].real(self.m_binsz)
                expcube["nxpix"].integer(self.m_nxpix)
                expcube["nypix"].integer(self.m_nypix)
                expcube["enumbins"].integer(self.m_ebins)
                expcube["emin"].real(self.m_emin)
                expcube["emax"].real(self.m_emax)   
                expcube["coordsys"].string(self.m_coordsys)
                expcube["proj"].string(self.m_proj)               
                expcube.run()
                
                # Header
                if self.logTerse():
                    self.log.header3("Creating PSF cube")
                
                # Create psf cube
                psfcube = ctools.ctpsfcube(select.obs())
                psfcube["incube"].filename("NONE")
                psfcube["usepnt"].boolean(False)
                psfcube["ebinalg"].string("LOG")
                psfcube["xref"].real(self.m_xref)
                psfcube["yref"].real(self.m_yref)
                psfcube["binsz"].real(self.m_binsz)
                psfcube["nxpix"].integer(self.m_nxpix)
                psfcube["nypix"].integer(self.m_nypix)
                psfcube["enumbins"].integer(self.m_ebins)
                psfcube["emin"].real(self.m_emin)
                psfcube["emax"].real(self.m_emax)    
                psfcube["coordsys"].string(self.m_coordsys)
                psfcube["proj"].string(self.m_proj)               
                psfcube.run()
                
                # Header
                if self.logTerse():
                    self.log.header3("Creating background cube")
                
                # Create background cube
                bkgcube = ctools.ctbkgcube(select.obs())
                bkgcube["incube"].filename("NONE")
                bkgcube["usepnt"].boolean(False)
                bkgcube["ebinalg"].string("LOG")
                bkgcube["xref"].real(self.m_xref)
                bkgcube["yref"].real(self.m_yref)
                bkgcube["binsz"].real(self.m_binsz)
                bkgcube["nxpix"].integer(self.m_nxpix)
                bkgcube["nypix"].integer(self.m_nypix)
                bkgcube["enumbins"].integer(self.m_ebins)
                bkgcube["emin"].real(self.m_emin)
                bkgcube["emax"].real(self.m_emax)   
                bkgcube["coordsys"].string(self.m_coordsys)
                bkgcube["proj"].string(self.m_proj)                
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
                    self.log("No event in this time bin. Bin is skipped\n")

                # Set all values to 0
                for col in columns:
                    col[i] = 0.0
                TSvalues[i]    = 0.0
                ulim_values[i] = 0.0
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
                ulimit["srcname"].string(self.m_srcname)
                ulimit["eref"].real(1.0)
                
                # Try to run upper limit and catch exceptions
                try:
                    ulimit.run()
                    ulimit_value = ulimit.flux_ulimit()
                except:
                    if self.logTerse():
                        self.log("Upper limit calculation failed\n")
                    ulimit_value = -1.0
            
            # Get TS value
            TS = -1.0
            if self.m_calc_ts:
                TS = source.ts() 
            
            # Set values for storage
            TSvalues[i] = TS
            
            # Set FITS column values
            for col in columns:
                if "e_" == col.name()[:2]:
                    col[i] = source.spectral()[col.name()[2:]].error()
                else:
                    col[i] = source.spectral()[col.name()].value()
            
            # Store upper limit value if available
            if ulimit_value > 0.0:
                ulim_values[i] = ulimit_value
         
            # Log information
            if self.logExplicit(): 
                self.log.header3("Results of bin "+str(i)+": MJD "+str(tmin.mjd())+"-"+str(tmax.mjd()))
                for col in columns:
                    if "e_" == col.name()[:2]:
                        continue
                    value = source.spectral()[col.name()].value()
                    error = source.spectral()[col.name()].error()
                    unit = source.spectral()[col.name()].unit()
                    self.log(" > "+col.name()+": "+str(value)+" +- "+str(error)+" "+unit+"\n")
                if self.m_calc_ts and TSvalues[i] > 0.0:
                    self.log(" > TS = "+str(TS)+" \n")
                if self.m_calc_ulimit and ulim_values[i] > 0.0:
                    self.log(" > UL = "+str(ulim_values[i])+" [ph/cm2/s]")
                self.log("\n")

        # Append filles columns to fits table    
        table.append(MJD)
        table.append(e_MJD)
        for col in columns:
            table.append(col)
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
    Generates light curve.
    """
    # Create instance of application
    app = cslightcrv(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
