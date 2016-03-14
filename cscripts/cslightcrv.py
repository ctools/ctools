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

    The cslightcrv class generates a light curve for Imaging Air Cherenkov
    Telescope event data by performing a maximum likelihood fit using
    ctlike in a series of time bins. The time bins can be either
    specified in an ASCII file, as an interval divided into equally
    sized time bins, or can be taken from the Good Time Intervals of the
    observation(s).

    The format of the ASCII file is one row per time bin, each specifying
    the start of stop value of the bin, separated by a whitespace. The
    times are given in Modified Julian Days (MJD). 
    
    Examples:
            >>> lcrv = cslightcrv()
            >>> lcrv.run()
            >>> ... (querying for parameters) ...
            >>> fits = lcrv.lightcurve()
                Generates a light curve and retrieves the results in
                a FITS file.

            >>> lcrv = cslightcrv()
            >>> lcrv.execute()
            >>> ... (querying for parameters) ...
                Generates a light curve and saves results in a FITS file.

            >>> lcrv = cslightcrv(obs)
            >>> lcrv.execute()
            >>> ... (querying for parameters) ...
                Generates a light curve from the observations in an
                observation container and saves results in a FITS file.
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
        self._srcname = ""
        self._tbins   = gammalib.GGti()
        self._ebins   = gammalib.GEbounds()
        self._binned  = False

        # Initialise observation
        if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
            self._obs = argv[0]
            argv      = argv[1:]
        else:      
            self._obs = gammalib.GObservations()

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
    def _get_parameters(self):
        """
        Get parameters from parfile.
        """
        # Set observation if not done before
        if self._obs == None or self._obs.size() == 0:
            self._require_inobs("cslightcrv::_get_parameters()")
            self._obs = self._get_observations()

        # Set models if we have none
        if self._obs.models().size() == 0:
            self._obs.models(self["inmodel"].filename())

        # Get source name   
        self._srcname = self["srcname"].string()

        # Get time boundaries             
        self._tbins = self._create_tbounds()

        # Unbinned or binned analysis?
        self._ebins    = self["enumbins"].integer()
        if self._ebins == 0:
            self._binned = False
        else:
            self._binned = True

        # Make sure that other parameters are queried now and not later
        self["emin"].real()
        self["emax"].real()
        if self._binned:
            self["coordsys"].string()
            self["proj"].string()
            self["xref"].real()
            self["yref"].real()
            self["nxpix"].integer()
            self["nypix"].integer()
            self["binsz"].real()
        
        # Do the same for the hidden parameters, just in case
        self["edisp"].boolean()
        self["calc_ulim"].boolean()
        self["calc_ts"].boolean()
        self["fix_bkg"].boolean()
        self["fix_srcs"].boolean()

        # Read ahead output parameters
        if self._read_ahead():
            self["outfile"].filename()

        #  Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Return
        return

    def _create_tbounds(self):
        """
        Creates light curve time bins.

        The method reads the following user parameters:
            tbinalg:  Time binning algorithm
            tbinfile: Time binning file (FITS or ASCII)
            tmin:     Start time (MJD)
            tmax:     Stop time (MJD)
            tbins:    Number of time bins

        Returns:
            Light curve bins in form of a GTI.
        """
        # Initialise Good Time Intervals
        gti = gammalib.GGti()

        # Get algorithm to use for defining the time intervals
        algorithm = self["tbinalg"].string()
        
        # Handle a FITS or a ASCII file for time bin definition
        if algorithm == "FILE":

            # Get the filename
            filename = self["tbinfile"].filename()
            
            # If the file a FITS file then load GTIs
            if filename.is_fits():
                gti.load(filename)

            # ... otherwise load file as CSV ASCII file
            csv = gammalib.GCsv(filename)
            for i in range(csv.nrows()):
                tmin = gammalib.GTime()
                tmax = gammalib.GTime()
                tmin.mjd(csv.real(i,0))
                tmax.mjd(csv.real(i,1))
                gti.append(tmin,tmax)

        # Handle linear time binning
        elif algorithm == "LIN":

            # Get start and stop time and number of time bins
            time_min = self["tmin"].real()
            time_max = self["tmax"].real()
            nbins    = self["tbins"].integer()            

            # Compute time step and setup time intervals
            time_step = (time_max - time_min) / float(nbins)
            for i in range(nbins):
                tmin = gammalib.GTime()
                tmax = gammalib.GTime()
                tmin.mjd(time_min +    i *time_step)
                tmax.mjd(time_min + (i+1)*time_step)
                gti.append(tmin,tmax)

        # Handle usage of observation GTIs
        elif algorithm == "GTI":

            # Append the GTIs of all observations
            for obs in self._obs:
                for i in range(obs.events().gti().size()):
                    gti.append(obs.events().gti().tstart(i),
                               obs.events().gti().tstop(i))
    
        # ... otherwise raise an exception (this should never occur)
        else:
            raise AttributeError('Paramter tbinalg="'+algorithm+'" unknown. '
                                 'Must be one of "FILE", "LIN" or "GTI".')

        # Return Good Time Intervals
        return gti

    def _get_free_par_names(self):
        """
        Return list of free parameter names.
        
        Returns:
            List of free parameter names.
        """
        # Initialise list of free parameter names
        names = []

        # Collect list of free parameter names
        for par in self._obs.models()[self._srcname]:
            if par.is_free():
                names.append(par.name())

        # Return names
        return names

    def _adjust_model_pars(self):
        """
        Adjust model parameters dependent on user parameters.
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Adjust model parameters")

        # Adjust model parameters dependent on input user parameters
        for model in self._obs.models():

            # Set TS flag for all models to false. The source of interest
            # will be set to true later
            model.tscalc(False)

            # Log model name
            if self._logNormal():
                self._log.header3(model.name())

            # Deal with the source of interest    
            if model.name() == self._srcname:
                if self["calc_ts"].boolean():
                    model.tscalc(True)

            elif (self["fix_bkg"].boolean() and
                  model.classname() != "GModelSky"):
                for par in model:
                    if par.is_free():
                        par.fix()
                        if self._logNormal():
                            self._log(gammalib.parformat(par.name()))
                            self._log("fixed\n")

            elif (self["fix_srcs"].boolean() and
                  model.classname() == "GModelSky"):
                for par in model:
                    if par.is_free():
                        par.fix()
                        if self._logNormal():
                            self._log(gammalib.parformat(par.name()))
                            self._log("fixed\n")

        # Return
        return

    def _create_fits_table(self, results):
        """
        Creates FITS binary table containing light curve results.

        Args:
            results: List of result dictionaries.

        Returns:
            FITS binary table containing light curve.
        """
        # Determine number of rows in FITS table
        nrows = len(results)

        # Create FITS Table with extension "LIGHTCURVE"
        table = gammalib.GFitsBinTable(nrows)
        table.extname("LIGHTCURVE")

        # Append time columns        
        mjd   = gammalib.GFitsTableDoubleCol("MJD", nrows)
        e_mjd = gammalib.GFitsTableDoubleCol("e_MJD", nrows)
        mjd.unit("days")
        e_mjd.unit("days")
        for i, result in enumerate(results):
            mjd[i]   = result["mjd"]
            e_mjd[i] = result["e_mjd"]
        table.append(mjd)
        table.append(e_mjd)

        # Create parameter columns
        for par in self._obs.models()[self._srcname]:
            if par.is_free():
                name   = par.name()
                e_name = "e_"+par.name()
                col    = gammalib.GFitsTableDoubleCol(name, nrows)
                e_col  = gammalib.GFitsTableDoubleCol(e_name, nrows)
                col.unit(par.unit())
                e_col.unit(par.unit())
                for i, result in enumerate(results):
                    col[i]   = result['values'][name]
                    e_col[i] = result['values'][e_name]
                table.append(col)
                table.append(e_col)

        # Append Test Statistic column
        ts = gammalib.GFitsTableDoubleCol("TS", nrows)
        for i, result in enumerate(results):
            ts[i] = result["ts"]
        table.append(ts)

        # Append upper limit column
        ulimit = gammalib.GFitsTableDoubleCol("UpperLimit", nrows)
        ulimit.unit("ph/cm2/s")
        for i, result in enumerate(results):
            ulimit[i] = result['ulimit']
        table.append(ulimit)

        # Return table
        return table

    def _compute_ulimit(self, obs):
        """
        Computes upper flux limit.

        Args:
            obs: Observation container.

        Returns:
            Upper flux limit (-1 of not computed).
        """
        # Initialise upper flux limit
        ulimit_value = -1.0
        
        # Perform computation only if requested
        if self["calc_ulim"].boolean():

            # Write header in logger
            if self._logExplicit():
                self._log.header3("Computing upper flux limit")

            # Create upper limit object  
            ulimit = ctools.ctulimit(obs)
            ulimit["srcname"] = self._srcname
            ulimit["eref"]    = 1.0

            # Try to run upper limit and catch exceptions
            try:
                ulimit.run()
                ulimit_value = ulimit.flux_ulimit()
            except:
                if self._logTerse():
                    self._log("Upper limit flux calculation failed.s\n")
                ulimit_value = -1.0

        # Return upper limit
        return ulimit_value

    def _bin_observation(self, obs):
        """
        Bin an observation is a binned analysis was requested.

        Args:
            obs: Observation container.

        Returns:
            Observation container with a binned.
        """
        # Header
        if self._logExplicit():
            self._log.header3("Binning events")

        # Bin events
        bin = ctools.ctbin(obs)
        bin["usepnt"]   = False
        bin["ebinalg"]  = "LOG"
        bin["xref"]     = self["xref"].real()
        bin["yref"]     = self["yref"].real()
        bin["binsz"]    = self["binsz"].real()
        bin["nxpix"]    = self["nxpix"].integer()
        bin["nypix"]    = self["nypix"].integer()
        bin["enumbins"] = self._ebins
        bin["emin"]     = self["emin"].real()
        bin["emax"]     = self["emax"].real()        
        bin["coordsys"] = self["coordsys"].string()
        bin["proj"]     = self["proj"].string()            
        bin.run()

        # Header
        if self._logExplicit():
            self._log.header3("Creating exposure cube")

        # Create exposure cube
        expcube = ctools.ctexpcube(obs)
        expcube["incube"]   = "NONE"
        expcube["usepnt"]   = False
        expcube["ebinalg"]  = "LOG"
        expcube["xref"]     = self["xref"].real()
        expcube["yref"]     = self["yref"].real()
        expcube["binsz"]    = self["binsz"].real()
        expcube["nxpix"]    = self["nxpix"].integer()
        expcube["nypix"]    = self["nypix"].integer()
        expcube["enumbins"] = self._ebins
        expcube["emin"]     = self["emin"].real()
        expcube["emax"]     = self["emax"].real()   
        expcube["coordsys"] = self["coordsys"].string()
        expcube["proj"]     = self["proj"].string()               
        expcube.run()

        # Header
        if self._logExplicit():
            self._log.header3("Creating point spread function cube")

        # Create psf cube
        psfcube = ctools.ctpsfcube(obs)
        psfcube["incube"]   = "NONE"
        psfcube["usepnt"]   = False
        psfcube["ebinalg"]  = "LOG"
        psfcube["xref"]     = self["xref"].real()
        psfcube["yref"]     = self["yref"].real()
        psfcube["binsz"]    = self["binsz"].real()
        psfcube["nxpix"]    = self["nxpix"].integer()
        psfcube["nypix"]    = self["nypix"].integer()
        psfcube["enumbins"] = self._ebins
        psfcube["emin"]     = self["emin"].real()
        psfcube["emax"]     = self["emax"].real()    
        psfcube["coordsys"] = self["coordsys"].string()
        psfcube["proj"]     = self["proj"].string()               
        psfcube.run()

        # Check if we need to include energy dispersion
        if self["edisp"].boolean():

            # Header
            if self._logExplicit():
                self._log.header3("Creating energy dispersion cube")
            
            # Create edisp cube
            edispcube = ctools.ctedispcube(obs)
            edispcube["incube"]   = "NONE"
            edispcube["usepnt"]   = False
            edispcube["ebinalg"]  = "LOG"
            edispcube["xref"]     = self["xref"].real()
            edispcube["yref"]     = self["yref"].real()
            edispcube["binsz"]    = self["binsz"].real()
            edispcube["nxpix"]    = self["nxpix"].integer()
            edispcube["nypix"]    = self["nypix"].integer()
            edispcube["enumbins"] = self._ebins
            edispcube["emin"]     = self["emin"].real()
            edispcube["emax"]     = self["emax"].real()    
            edispcube["coordsys"] = self["coordsys"].string()
            edispcube["proj"]     = self["proj"].string()               
            edispcube.run()

        # Header
        if self._logExplicit():
            self._log.header3("Creating background cube")

        # Create background cube
        bkgcube = ctools.ctbkgcube(obs)
        bkgcube["incube"]   = "NONE"
        bkgcube["usepnt"]   = False
        bkgcube["ebinalg"]  = "LOG"
        bkgcube["xref"]     = self["xref"].real()
        bkgcube["yref"]     = self["yref"].real()
        bkgcube["binsz"]    = self["binsz"].real()
        bkgcube["nxpix"]    = self["nxpix"].integer()
        bkgcube["nypix"]    = self["nypix"].integer()
        bkgcube["enumbins"] = self._ebins
        bkgcube["emin"]     = self["emin"].real()
        bkgcube["emax"]     = self["emax"].real()   
        bkgcube["coordsys"] = self["coordsys"].string()
        bkgcube["proj"]     = self["proj"].string()                
        bkgcube.run()

        # Retrieve a new oberservation container
        new_obs = bin.obs().copy()
        
        # Get new models
        models = bkgcube.models()
        
        # Set stacked response
        if self["edisp"].boolean():
            new_obs[0].response(expcube.expcube(),
                                psfcube.psfcube(),
                                edispcube.edispcube(),
                                bkgcube.bkgcube())                    
        else:
            new_obs[0].response(expcube.expcube(),
                                psfcube.psfcube(),
                                bkgcube.bkgcube())

        # Fix background models if required
        if self["fix_bkg"].boolean():
            for model in models:
                if model.classname() != "GModelSky":
                    for par in model:
                        par.fix()

        # Set models for new oberservation container     
        new_obs.models(models)

        # Return new oberservation container
        return new_obs


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

        # Adjust model parameters dependent on user parameters
        self._adjust_model_pars()

        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Generate lightcurve")      

        # Initialise list of result dictionaries
        results = []

        # Get source parameters
        pars = self._get_free_par_names()

        # Loop over time bins
        for i in range(self._tbins.size()):

            # Get time boundaries
            tmin = self._tbins.tstart(i)
            tmax = self._tbins.tstop(i)

            # Write time bin into header
            if self._logTerse():
                self._log.header2("MJD "+
                                  str(tmin.mjd())+"-"+
                                  str(tmax.mjd()))

            # Compute time bin center and time width
            twidth = 0.5 * (tmax - tmin) # in seconds
            tmean  = tmin + twidth

            # Initialise result dictionary
            result = {'mjd': tmean.mjd(),
                      'e_mjd': twidth / gammalib.sec_in_day,
                      'ts': 0.0,
                      'ulimit': 0.0,
                      'pars': pars,
                      'values': {}}

            # Log information
            if self._logExplicit():
                self._log.header3("Selecting events")

            # Select events
            select = ctools.ctselect(self._obs)
            select["emin"] = self["emin"].real()    
            select["emax"] = self["emax"].real() 
            select["tmin"] = tmin.convert(select._time_reference())
            select["tmax"] = tmax.convert(select._time_reference())
            select["rad"]  = "UNDEFINED"
            select["ra"]   = "UNDEFINED"
            select["dec"]  = "UNDEFINED"
            select.run()  

            # Retrieve observation
            obs = select.obs()

            # If a stacked analysis is requested then bin the events
            # and compute the stacked response functions and setup
            # an observation container with a single stacked observation.
            if self._binned:
                obs = self._bin_observation(obs)

            # Header
            if self._logExplicit():
                self._log.header3("Fitting the data")

            # Do maximum likelihood model fitting
            like = ctools.ctlike(obs)
            like["edisp"] = self["edisp"].boolean()
            like.run()

            # Skip bin if no event was present
            if like.obs().logL() == 0.0:

                # Signal skipping of bin
                if self._logTerse():
                    self._log(gammalib.parformat("Warning"))
                    self._log("No event in this time bin, skip bin.\n")

                # Set all results to 0
                for par in pars:
                    result['values'][par]      = 0.0
                    result['values']["e_"+par] = 0.0

                # Append result
                results.append(result)

                # Continue with next time bin
                continue

            # Retrieve model fitting results for source of interest
            source = like.obs().models()[self._srcname]

            # Extract parameter values
            for par in pars:
                result['values'][par]      = source.spectral()[par].value()
                result['values']["e_"+par] = source.spectral()[par].error()

            # Calculate upper limit (-1 if not computed)
            ulimit_value = self._compute_ulimit(like.obs())
            if ulimit_value > 0.0:
                result['ulimit'] = ulimit_value

            # Extract Test Statistic value
            if self["calc_ts"].boolean():
                result['ts'] = source.ts() 

            # Append result to list of dictionaries
            results.append(result)

            # Log results for this time bin
            if self._logNormal():
                self._log.header3("Results")
                pars = self._get_free_par_names()
                for par in pars:
                    value = source.spectral()[par].value()
                    error = source.spectral()[par].error()
                    unit  = source.spectral()[par].unit()
                    self._log(gammalib.parformat(par))
                    self._log(str(value))
                    self._log(" +/- ")
                    self._log(str(error))
                    self._log(" ")
                    self._log(unit)
                    self._log("\n")
                if result['ulimit'] > 0.0:
                    self._log(gammalib.parformat("Upper flux limit"))
                    self._log(str(result['ulimit'])+" ph/cm2/s\n")
                if self["calc_ts"].boolean():
                    self._log(gammalib.parformat("Test Statistic"))
                    self._log(str(result['ts'])+"\n")

        # Create FITS table from results
        table = self._create_fits_table(results)

        # Create FITS file and append FITS table
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

        # Save lightcurve
        self.save()

        # Return
        return

    def save(self):
        """
        Save light curve.
        """
        # Write header
        if self._logTerse():
            self._log('\n')
            self._log.header1('Save light curve')

        # Get light curve filename
        outfile = self["outfile"].filename()
        
        # Continue only filename and residual map are valid
        if self._fits != None:

            # Log file name
            if self._logTerse():
                self._log(gammalib.parformat("Light curve file"))
                self._log(outfile.url())
                self._log("\n")

            # Save spectrum
            self._fits.saveto(outfile, self._clobber())

        # Return
        return

    def lightcurve(self):
        """
        Return light curve FITS file.

        Returns:
            FITS file containing light curve.
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

    # Execute application
    app.execute()
