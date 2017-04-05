#! /usr/bin/env python
# ==========================================================================
# Generates pulsar spectra as a function of phase
#
# Copyright (C) 2014-2017 Rolf Buehler
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
import sys
import gammalib
import ctools
from cscripts import obsutils


# ================ #
# csphasecrv class #
# ================ #
class csphasecrv(ctools.cscript):
    """
    Generates a spectra in phase bins

    This script computes spectra by performing a maximum likelihood fit
    using :doc:`ctlike` in a series of phase bins for pulsars.
    The phase bins can be either specified in an ASCII file, as an interval
    divided into equally sized phase bins. 

    Examples:
            >>> phrv = csphasecrv()
            >>> phrv.run()
            >>> ... (querying for parameters) ...
            >>> phrv = phrv.phasecrv()
                Generates phase fits and retrieves dictionary with best fit models.

            >>> lcrv = csphasecrv(obs)
            >>> lcrv.execute()
            >>> ... (querying for parameters) ...
                Generates phase fits from the observations 
                and saves results in XML files.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'csphasecrv'
        self._version = '1.0.0'

        # Initialise some members
        self._srcname = ''
        # Phases are stored in a nested list [[ph1min,ph1max], [ph2min,ph2max],..]
        self._phbins  = [[0.0,1.0]]   
        self._stacked = False

        # Initialise observation container from constructor arguments.
        self._obs, argv = self._set_input_obs(argv)

        # Initialise script by calling the appropriate class constructor.
        self._init_cscript(argv)

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Setup observations (require response and allow event list, don't
        # allow counts cube)
        self._setup_observations(self._obs, True, True, False)

        # Set models if there are none in the container
        if self._obs.models().size() == 0:
            self._obs.models(self['inmodel'].filename())

        # Get source name
        self._srcname = self['srcname'].string()

        # Get time boundaries
        self._phbins = self._create_tbounds()

        # Set stacked analysis flag to True if the requested number of
        # energy bins is positive. Otherwise an unbinned analysis will
        # be done and the stacked analysis flag will be set to False.
        if self['enumbins'].integer() > 0:
            self._stacked = True
        else:
            self._stacked = False

        # Make sure that remaining user parameters are queried now. We
        # do not store the actual parameter values as we do not want
        # too many instance attributes with enhances the maintenance
        # costs.
        self['emin'].real()
        self['emax'].real()
        if self._stacked:
            self['coordsys'].string()
            self['proj'].string()
            self['xref'].real()
            self['yref'].real()
            self['nxpix'].integer()
            self['nypix'].integer()
            self['binsz'].real()

        # Do the same for the hidden parameters, just in case
        self['edisp'].boolean()
        
        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _create_tbounds(self):
        """
        Creates phase bins

        Returns
        -------
        phbins : `[]`
            List of phase bins
        """
        # Initialise Good Time Intervals
        phbins = []

        # Get algorithm to use for defining the time intervals. Possible
        # values are "FILE" or "LIN". This is enforced at
        # parameter file level, hence no checking is needed.
        algorithm = self['tbinalg'].string()

        # If the algorithm is "FILE" then handle a FITS or an ASCII file for
        # the time bin definition
        if algorithm == 'FILE':

            # Get the filename
            filename = self['phbinfile'].filename()

            # Load  ASCII file as CSV file and construct
            # the GTIs from the rows of the CSV file. It is expected that the
            # CSV file has two columns containing the "START" and "STOP"
            # values of the phase bins. No header row is expected.
            csv = gammalib.GCsv(filename)
            for i in range(csv.nrows()):
                phmin = csv.real(i,0)
                phmax = csv.real(i,1)
                phbins.append([phmin,phmax])

        # If the algorithm is "LIN" then use a linear time binning, defined by
        # the "phbins" user parameters
        elif algorithm == 'LIN':

            # Get start and stop time and number of time bins
            nbins    = self['phbins'].integer()

            # Compute time step and setup the GTIs
            ph_step = 1.0 / float(nbins)
            for i in range(nbins):
                phmin = i *ph_step
                phmax = (i+1)*ph_step
                phbins.append([phmin,phmax])

        # Return Good Time Intervals
        return phbins

    def _bin_observation(self, obs):
        """
        Bin an observation if a binned analysis was requested

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Observation container

        Returns
        -------
        obs : `~gammalib.GObservations`
            Observation container where the first observation is a binned observation
        """
        # Header
        if self._logExplicit():
            self._log.header3('Binning events')

        # Bin events
        cntcube = ctools.ctbin(obs)
        cntcube['usepnt']   = False
        cntcube['ebinalg']  = 'LOG'
        cntcube['xref']     = self['xref'].real()
        cntcube['yref']     = self['yref'].real()
        cntcube['binsz']    = self['binsz'].real()
        cntcube['nxpix']    = self['nxpix'].integer()
        cntcube['nypix']    = self['nypix'].integer()
        cntcube['enumbins'] = self['enumbins'].integer()
        cntcube['emin']     = self['emin'].real()
        cntcube['emax']     = self['emax'].real()
        cntcube['coordsys'] = self['coordsys'].string()
        cntcube['proj']     = self['proj'].string()
        cntcube.run()

        # Header
        if self._logExplicit():
            self._log.header3('Creating stacked response')

        # Get stacked response
        response = obsutils.get_stacked_response(obs,
                                                 self['xref'].real(),
                                                 self['yref'].real(),
                                                 binsz=self['binsz'].real(),
                                                 nxpix=self['nxpix'].integer(),
                                                 nypix=self['nypix'].integer(),
                                                 emin=self['emin'].real(),
                                                 emax=self['emax'].real(),
                                                 enumbins=self['enumbins'].integer(),
                                                 edisp=self['edisp'].boolean(),
                                                 coordsys=self['coordsys'].string(),
                                                 proj=self['proj'].string())

        # Retrieve a new oberservation container
        new_obs = cntcube.obs().copy()

        # Set stacked response
        if self['edisp'].boolean():
            new_obs[0].response(response['expcube'], response['psfcube'],
                                response['edispcube'], response['bkgcube'])
        else:
            new_obs[0].response(response['expcube'], response['psfcube'],
                                response['bkgcube'])

        # Get new models
        models = response['models']

        # Set models for new oberservation container     
        new_obs.models(models)

        # Return new oberservation container
        return new_obs


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)
        
        # Get parameters
        self._get_parameters()
        
        # Write observation into logger
        if self._logTerse():
            self._log('\n')
            self._log.header1(gammalib.number('Input observation',len(self._obs)))
            self._log(str(self._obs))
            self._log('\n')
                
        #Store and clear observations to split them in phase
        orig_obs = self._obs.copy()
        #self._obs.clear()
        
        #Dictionary to save phase fitted models
        self._fitmodels = {}
        
        #Loop over all phases
        for phbin in self._phbins:
            
            # Select events
            select = ctools.ctselect(orig_obs)
            select['emin'] = self['emin'].real()
            select['emax'] = self['emax'].real()
            select['tmin'] = 'UNDEFINED'
            select['tmax'] = 'UNDEFINED'
            select['rad']  = 'UNDEFINED'
            select['ra']   = 'UNDEFINED'
            select['dec']  = 'UNDEFINED'
            select['expr'] = "PHASE>"+str(phbin[0])+" && PHASE<"+str(phbin[1])
            select.run()
            
            #Add phase to observation id
            for i in range(0,select.obs().size()):
                oldid = select.obs()[i].id()
                select.obs()[i].id(oldid+"_"+str(phbin[0])+"-"+str(phbin[1]))
            obs = select.obs()
            
            # If a stacked analysis is requested then bin the events
            # and compute the stacked response functions and setup
            # an observation container with a single stacked observation.
            if self._stacked:
                obs = self._bin_observation(select.obs())
    
            # Header
            if self._logExplicit():
                self._log.header3('Fitting the data')
    
            # Do maximum likelihood model fitting
            like = ctools.ctlike(obs)
            like['edisp'] = self['edisp'].boolean()
            like.run()
            
            # Renormalize models to phase selection
            # TODO move the scaling from the temporal to the spectral component
            for model in like.obs().models():
                scaled_norm = model["Normalization"].value()/(phbin[1] - phbin[0])
                model["Normalization"].value(scaled_norm)
            
            #Store fit model
            self._fitmodels[str(phbin[0])+"-"+str(phbin[1])] = like.obs().models().copy()

            ## Store observations
            #for ob in select.obs():
                #self._obs.append(ob.copy())

        # Return
        return

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save fitted models
        self.save()

        # Return
        return

    def saveXML(self):
        """
        Save phase fits into one XML model each
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save phase XML')

        # Get out filename
        outfile = str()

        # Save XMLs
        for phbin in self._phbins:
            phname = str(phbin[0])+"-"+str(phbin[1])
            outname = self['outfile'].filename()
            outname_ph = (outname.file()).replace(".fits","_ph_"+phname+".xml")
            outfile = outname.path()+outname_ph
            print outfile
            self._fitmodels[phname].save(outfile)
        
        # Return
        return
    
    def save(self):
        """
        Saves results to fits and XML
        """
        print "SAVING"
        self.savefits()
        self.saveXML()
        
        # Return
        return

    def phaserv(self):
        """
        Return dictionary with best fit models
        """
        # Return
        return self._fitmodels

    def _get_free_par_names(self):
        """
        Return list of free parameter names

        Returns
        -------
        names : list of str
            List of free parameter names.
        """
        # Initialise empty list of free parameter names
        names = []

        # Collect list of free parameter names
        for par in self._obs.models()[self._srcname]:
            if par.is_free():
                names.append(par.name())

        # Return names
        return names

    def savefits(self):
        """
        Saved phase fit results to a fits file
        """
        # Initialise list of result dictionaries
        results = []

        # Get source parameters
        pars = self._get_free_par_names()

        # Loop over time bins
        for i in range(len(self._phbins)):

            # Get time boundaries
            phmin = self._phbins[i][0]
            phmax = self._phbins[i][1]

            # Write time bin into header
            if self._logTerse():
                self._log.header2('PHASE '+str(phmin)+'-'+str(phmax))

            # Initialise result dictionary
            result = {'phmin': phmin,
                      'phmax': phmax,
                      'pars': pars,
                      'values': {}}
                        # Extract parameter values
            
            #Store fit results
            phname = str(phmin)+"-"+str(phmax)
            source = self._fitmodels[phname][self._srcname]
            for par in pars:
                result['values'][par]      = source[par].value()
                result['values']['e_'+par] = source[par].error()
            
            # Append result to list of dictionaries
            results.append(result)
            
        # Create FITS table from results
        table = self._create_fits_table(results)

        # Create FITS file and append FITS table to FITS file
        self._fits = gammalib.GFits()
        self._fits.append(table)

        # Save to fits
        outfile = self['outfile'].filename()
        print outfile
        self._fits.saveto(outfile, self._clobber())
        
        #Return
        return

    def _create_fits_table(self, results):
        """
        Creates FITS binary table containing light curve results

        Parameters
        ----------
        results : list of dict
            List of result dictionaries

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            FITS binary table containing light curve
        """
        # Determine number of rows in FITS table
        nrows = len(results)

        # Create FITS Table with extension "LIGHTCURVE"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('PHASESPEC')

        # Append time columns "MJD" and "e_MJD"
        phmin   = gammalib.GFitsTableDoubleCol('PHMIN', nrows)
        phmax = gammalib.GFitsTableDoubleCol('PHMAX', nrows)
        #mjd.unit('days')
        #e_mjd.unit('days')
        for i, result in enumerate(results):
            phmin[i]   = result['phmin']
            phmax[i]   = result['phmax']
        table.append(phmin)
        table.append(phmax)

        # Create parameter columns
        for par in self._obs.models()[self._srcname]:
            if par.is_free():
                name   = par.name()
                e_name = 'e_'+par.name()
                col    = gammalib.GFitsTableDoubleCol(name, nrows)
                e_col  = gammalib.GFitsTableDoubleCol(e_name, nrows)
                col.unit(par.unit())
                e_col.unit(par.unit())
                for i, result in enumerate(results):
                    col[i]   = result['values'][name]
                    e_col[i] = result['values'][e_name]
                table.append(col)
                table.append(e_col)

        # Return table
        return table

        
        
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csphasecrv(sys.argv)

    # Execute application
    app.execute()
