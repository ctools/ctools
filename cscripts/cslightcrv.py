#! /usr/bin/env python
# ==========================================================================
# Generates a lightcurve.
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


# ================ #
# cslightcrv class #
# ================ #
class cslightcrv(ctools.cscript):
    """
    Generates a lightcurve.

    This class implements the creation of a light curve. It derives from
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
        self._name    = "cslightcrv"
        self._version = "1.1.0"

        # Initialise some members
        self._outfile = ""

        # Initialise observation
        if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
            self._obs = argv[0]
            argv      = argv[1:]
        else:      
            self._obs = gammalib.GObservations()

        # Make sure that parfile exists
        self._parfile()

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self._name, self._version)
        elif len(argv) ==1:
            ctools.cscript.__init__(self, self._name, self._version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self._log_header()
        self._log.date(True)

        # Return
        return

    def __del__(self):
        """
        Destructor.
        """
        # Return
        return


    # Private methods
    def _parfile(self):
        """
        Check if parfile exists. If parfile does not exist then create a
        default parfile. This kluge avoids shipping the cscript with a parfile.
        """
        # Set parfile name
        parfile = self._name+".par"

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
            pars.append(gammalib.GApplicationPar("edispcube","f","a","NONE","","","Input energy dispersion cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("bkgcube","s","a","NONE","","","Input background cube file (only needed for stacked analysis)"))
            pars.append(gammalib.GApplicationPar("edisp","b","h","no","","","Apply energy dispersion?"))
            pars.append(gammalib.GApplicationPar("caldb","s","a","prod2","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","South_0.5h","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","lightcurve.fits","","","Output light curve file"))
            pars.append(gammalib.GApplicationPar("tbinalg","s","a","GTI","FILE|LIN|GTI","", "Algorithm for defining time bins"))
            pars.append(gammalib.GApplicationPar("tmin","r","a","51544.5","","days","Lightcurve start time [MJD]"))
            pars.append(gammalib.GApplicationPar("tmax","r","a","51544.6","","days","Lightcurve stop time [MJD]"))
            pars.append(gammalib.GApplicationPar("tbins","i","a","5","","","Number of time bins"))
            pars.append(gammalib.GApplicationPar("tbinfile","f","a","tbins.fits","","", "File defining the time binning"))
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

    def _create_tbounds(self):

        # Create time bin container
        self._tbins = gammalib.GGti()

        # check for binning algorithm
        if self["tbinalg"].string() == "FILE":
            try:
                self._tbins.load(self["tbinfile"].filename())
            except:               
                csv = gammalib.GCsv(self["tbinfile"].filename())
                for i in range(csv.nrows()):
                    tmin = gammalib.GTime()
                    tmax = gammalib.GTime()
                    tmin.mjd(csv.real(i,0))
                    tmax.mjd(csv.real(i,1))
                    self._tbins.append(tmin,tmax)
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
                self._tbins.append(tmin,tmax)
        elif self["tbinalg"].string() == "GTI":
            # Use GTIs of observations
            for obs in self._obs:
                for i in range(obs.events().gti().size()):
                    self._tbins.append(obs.events().gti().tstart(i),obs.events().gti().tstop(i))
        else:
            raise AttributeError("tbinalg=\""+self["tbinalg"].string()+"\" unkown. Must be one of \"FILE\", \"LIN\" or \"GTI\"")

        # Return
        return

    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Set observation if not done before
        if self._obs == None or self._obs.size() == 0:
            self._require_inobs("cslightcrv::get_parameters()")
            self._obs = self._get_observations()

        # Set models if we have none
        if self._obs.models().size() == 0:
            self._obs.models(self["inmodel"].filename())

        # Get source name   
        self._srcname = self["srcname"].string()

        # Get energy dispersion flag
        self._edisp = self["edisp"].boolean()
        
        # Get time boundaries             
        self._create_tbounds()

        # Unbinned or binned analysis?
        self._ebins    = self["enumbins"].integer()
        if self._ebins == 0:
            self._binned = False
        else:
            self._binned = True

        # Get energy range
        self._emin = self["emin"].real()
        self._emax = self["emax"].real()

        # Get binning flag
        if self._binned:
            self._coordsys = self["coordsys"].string()
            self._proj     = self["proj"].string()
            self._xref     = self["xref"].real()
            self._yref     = self["yref"].real()
            self._nxpix    = self["nxpix"].integer()
            self._nypix    = self["nypix"].integer()
            self._binsz    = self["binsz"].real()

        # Read other parameters
        self._outfile = self["outfile"].filename()

        # Get other parameeters
        self._calc_ulimit = self["calc_ulim"].boolean()
        self._calc_ts     = self["calc_ts"].boolean()
        self._fix_bkg     = self["fix_bkg"].boolean()
        self._fix_srcs    = self["fix_srcs"].boolean()

        # Set some fixed parameters
        self._chatter = self["chatter"].integer()
        self._clobber = self["clobber"].boolean()
        self._debug   = self["debug"].boolean()

        #  Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

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
            self._log.header1("Generate lightcurve")      

        # Initialise FITS Table with extension "LIGHTCURVE"
        table = gammalib.GFitsBinTable(self._tbins.size())
        table.extname("LIGHTCURVE")

        # Add Header for compatibility with gammalib.GMWLSpectrum
        table.card("INSTRUME", "CTA", "Name of Instrument")
        table.card("TELESCOP", "CTA", "Name of Telescope")

        # Create FITS table columns        
        MJD = gammalib.GFitsTableDoubleCol("MJD", self._tbins.size())
        MJD.unit("days")
        e_MJD = gammalib.GFitsTableDoubleCol("e_MJD", self._tbins.size())
        e_MJD.unit("days")

        # Create a FITS column for every free parameter
        columns = []
        for par in self._obs.models()[self._srcname]:
            if par.is_free():
                col = gammalib.GFitsTableDoubleCol(par.name(), self._tbins.size())
                col.unit(par.unit())
                columns.append(col)
                e_col = gammalib.GFitsTableDoubleCol("e_"+par.name(), self._tbins.size())
                e_col.unit(par.unit())
                columns.append(e_col)

        # Create TS and upper limit columns
        TSvalues    = gammalib.GFitsTableDoubleCol("TS", self._tbins.size())
        ulim_values = gammalib.GFitsTableDoubleCol("UpperLimit", self._tbins.size())
        ulim_values.unit("ph/cm2/s")

        # Loop over energy bins
        for i in range(self._tbins.size()):

            # Log information
            if self._logTerse():
                self._log("\n")
                self._log.header2("Time bin "+str(i))

            # Get time boundaries
            tmin = self._tbins.tstart(i)
            tmax = self._tbins.tstop(i)

            # Compute time bin center and time width
            twidth = 0.5 * (tmax - tmin) # in seconds
            tmean  = tmin + twidth

            # Store time as MJD
            MJD[i]   = tmean.mjd()
            e_MJD[i] = twidth / gammalib.sec_in_day # in days

            # Log information
            if self._logExplicit():
                self._log.header3("Selecting events")

            # Select events
            select = ctools.ctselect(self._obs)
            select["emin"] = self._emin    
            select["emax"] = self._emax 
            select["tmin"] = tmin.convert(select._time_reference())
            select["tmax"] = tmax.convert(select._time_reference())
            select["rad"]  = "UNDEFINED"
            select["ra"]   = "UNDEFINED"
            select["dec"]  = "UNDEFINED"
            select.run()  

            # Retrieve observation
            obs = select.obs()

            # Binned analysis
            if self._binned:

                # Header
                if self._logTerse():
                    self._log.header3("Binning events")

                # Bin events
                bin = ctools.ctbin(select.obs())
                bin["usepnt"]   = False
                bin["ebinalg"]  = "LOG"
                bin["xref"]     = self._xref
                bin["yref"]     = self._yref
                bin["binsz"]    = self._binsz
                bin["nxpix"]    = self._nxpix
                bin["nypix"]    = self._nypix
                bin["enumbins"] = self._ebins
                bin["emin"]     = self._emin
                bin["emax"]     = self._emax        
                bin["coordsys"] = self._coordsys
                bin["proj"]     = self._proj            
                bin.run()

                # Header
                if self._logTerse():
                    self._log.header3("Creating exposure cube")

                # Create exposure cube
                expcube = ctools.ctexpcube(select.obs())
                expcube["incube"]   = "NONE"
                expcube["usepnt"]   = False
                expcube["ebinalg"]  = "LOG"
                expcube["xref"]     = self._xref
                expcube["yref"]     = self._yref
                expcube["binsz"]    = self._binsz
                expcube["nxpix"]    = self._nxpix
                expcube["nypix"]    = self._nypix
                expcube["enumbins"] = self._ebins
                expcube["emin"]     = self._emin
                expcube["emax"]     = self._emax   
                expcube["coordsys"] = self._coordsys
                expcube["proj"]     = self._proj               
                expcube.run()

                # Header
                if self._logTerse():
                    self._log.header3("Creating point spread function cube")

                # Create psf cube
                psfcube = ctools.ctpsfcube(select.obs())
                psfcube["incube"]   = "NONE"
                psfcube["usepnt"]   = False
                psfcube["ebinalg"]  = "LOG"
                psfcube["xref"]     = self._xref
                psfcube["yref"]     = self._yref
                psfcube["binsz"]    = self._binsz
                psfcube["nxpix"]    = self._nxpix
                psfcube["nypix"]    = self._nypix
                psfcube["enumbins"] = self._ebins
                psfcube["emin"]     = self._emin
                psfcube["emax"]     = self._emax    
                psfcube["coordsys"] = self._coordsys
                psfcube["proj"]     = self._proj               
                psfcube.run()

                # Check if we need to include energy dispersion
                if self._edisp:

                    # Header
                    if self._logTerse():
                        self._log.header3("Creating energy dispersion cube")
                    
                    # Create edisp cube
                    edispcube = ctools.ctedispcube(select.obs())
                    edispcube["incube"]   = "NONE"
                    edispcube["usepnt"]   = False
                    edispcube["ebinalg"]  = "LOG"
                    edispcube["xref"]     = self._xref
                    edispcube["yref"]     = self._yref
                    edispcube["binsz"]    = self._binsz
                    edispcube["nxpix"]    = self._nxpix
                    edispcube["nypix"]    = self._nypix
                    edispcube["enumbins"] = self._ebins
                    edispcube["emin"]     = self._emin
                    edispcube["emax"]     = self._emax    
                    edispcube["coordsys"] = self._coordsys
                    edispcube["proj"]     = self._proj               
                    edispcube.run()

                # Header
                if self._logTerse():
                    self._log.header3("Creating background cube")

                # Create background cube
                bkgcube = ctools.ctbkgcube(select.obs())
                bkgcube["incube"]   = "NONE"
                bkgcube["usepnt"]   = False
                bkgcube["ebinalg"]  = "LOG"
                bkgcube["xref"]     = self._xref
                bkgcube["yref"]     = self._yref
                bkgcube["binsz"]    = self._binsz
                bkgcube["nxpix"]    = self._nxpix
                bkgcube["nypix"]    = self._nypix
                bkgcube["enumbins"] = self._ebins
                bkgcube["emin"]     = self._emin
                bkgcube["emax"]     = self._emax   
                bkgcube["coordsys"] = self._coordsys
                bkgcube["proj"]     = self._proj                
                bkgcube.run()

                # Set new binned observation
                obs = bin.obs()
                
                # Get new models
                models = bkgcube.models()
                
                # Set precomputed binned response
                if self._edisp:
                    obs[0].response(expcube.expcube(), psfcube.psfcube(),
                                    edispcube.edispcube(), bkgcube.bkgcube())                    
                else:
                    obs[0].response(expcube.expcube(), psfcube.psfcube(),
                                    bkgcube.bkgcube())

                # Fix background models if required
                if self._fix_bkg:
                    for model in models:
                        if not model.classname() == "GModelSky":
                            for par in model:
                                par.fix()

                # Set new models to binned observation           
                obs.models(models)

            # Header
            if self._logTerse():
                self._log.header3("Performing fit")

            # Likelihood
            like = ctools.ctlike(obs)
            like["edisp"] = self._edisp
            like.run()

            # Skip bin if no event was present
            if like.obs().logL() == 0.0:

                # Log information
                if self._logTerse():
                    self._log("No event in this time bin. Bin is skipped\n")

                # Set all values to 0
                for col in columns:
                    col[i] = 0.0
                TSvalues[i]    = 0.0
                ulim_values[i] = 0.0
                continue

            # Get results
            fitted_models = like.obs().models()
            source        = fitted_models[self._srcname]

            # Calculate Upper Limit            
            ulimit_value = -1.0
            if self._calc_ulimit:

                # Logging information
                if self._logTerse():
                    self._log.header3("Computing upper limit")

                # Create upper limit object  
                ulimit = ctools.ctulimit(like.obs())
                ulimit["srcname"] = self._srcname
                ulimit["eref"] = 1.0

                # Try to run upper limit and catch exceptions
                try:
                    ulimit.run()
                    ulimit_value = ulimit.flux_ulimit()
                except:
                    if self._logTerse():
                        self._log("Upper limit calculation failed\n")
                    ulimit_value = -1.0

            # Get TS value
            TS = -1.0
            if self._calc_ts:
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
            if self._logExplicit(): 
                self._log.header3("Results of bin "+str(i)+": MJD "+str(tmin.mjd())+"-"+str(tmax.mjd()))
                for col in columns:
                    if "e_" == col.name()[:2]:
                        continue
                    value = source.spectral()[col.name()].value()
                    error = source.spectral()[col.name()].error()
                    unit = source.spectral()[col.name()].unit()
                    self._log(" > "+col.name()+": "+str(value)+" +- "+str(error)+" "+unit+"\n")
                if self._calc_ts and TSvalues[i] > 0.0:
                    self._log(" > TS = "+str(TS)+" \n")
                if self._calc_ulimit and ulim_values[i] > 0.0:
                    self._log(" > UL = "+str(ulim_values[i])+" [ph/cm2/s]")
                self._log("\n")

        # Append filles columns to fits table    
        table.append(MJD)
        table.append(e_MJD)
        for col in columns:
            table.append(col)
        table.append(TSvalues)
        table.append(ulim_values)

        # Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append(table)

        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save lightcurve
        self.save()

        # Return
        return

    def save(self):
        """
        Save lightcurve.
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save lightcurve")

        # Get outmap parameter
        self._outfile = self["outfile"].filename()
        
        # Continue only filename and residual map are valid
        if self._outfile != "NONE" and self._fits != None:

            # Log file name
            if self._logTerse():
                self._log("Save lightcurve into \""+self._outfile+"\".\n")

            # Save spectrum
            self._fits.saveto(self._outfile, self._clobber)

        # Return
        return

    def lightcurve(self):
        """
        Return lightcurve FITS file.

        Returns:
            FITS file containing lightcurve.
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
    app = cslightcrv(sys.argv)

    # Open logfile
    app._logFileOpen()

    # Execute application
    app.execute()
