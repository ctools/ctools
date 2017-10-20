#! /usr/bin/env python
# ==========================================================================
# Generates spectra as a function of phase
#
# Copyright (C) 2017 Rolf Buehler
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
from cscripts import ioutils


# ================ #
# csphasecrv class #
# ================ #
class csphasecrv(ctools.cscript):
    """
    Generates spectra in phase bins

    This script computes spectra by performing a maximum likelihood fit
    using :doc:`ctlike` in a series of phase bins for pulsars. The phase bins
    can be either specified in an ASCII file or as an interval divided into
    equally sized phase bins.

    Examples:
            >>> phcrv = csphasecrv()
            >>> phcrv.run()
            >>> ... (querying for parameters) ...
            >>> phcrv = phrv.phasecrv()
                Generates phase fits and retrieves dictionary with best fit models.

            >>> phcrv = csphasecrv(obs)
            >>> phcrv.execute()
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
        self._version = ctools.__version__

        # Initialise some members. Phases are stored in a nested list
        # [[ph1min,ph1max], [ph2min,ph2max],..]
        self._srcname = ''
        self._phbins  = [[0.0,1.0]]
        self._stacked = False
        self._fits    = gammalib.GFits()

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

        # Get phase boundaries
        self._phbins = self._create_tbounds()

        # Set stacked analysis flag and query relevant user parameters
        self._stacked = self._is_stacked()

        # Query the hidden parameters, just in case
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
        algorithm = self['phbinalg'].string()

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
                phmin = i     * ph_step
                phmax = (i+1) * ph_step
                phbins.append([phmin,phmax])

        # Return Good Time Intervals
        return phbins

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

    def _create_fits_table(self, results):
        """
        Creates FITS binary table containing phase curve results

        Parameters
        ----------
        results : list of dict
            List of result dictionaries

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            FITS binary table containing phase curve
        """
        # Determine number of rows in FITS table
        nrows = len(results)

        # Create FITS Table with extension "PHASECURVE"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('PHASECURVE')

        # Append phase columns
        ioutils.append_result_column(table, results, 'PHASE_MIN', '', 'phmin')
        ioutils.append_result_column(table, results, 'PHASE_MAX', '', 'phmax')

        # Append parameter columns
        ioutils.append_model_par_column(table, self._obs.models()[self._srcname],
                                        results)

        # Return table
        return table

    def _create_fits(self):
        """
        Create FITS file object from fit results
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

            # Initialise result dictionary
            result = {'phmin': phmin,
                      'phmax': phmax,
                      'pars': pars,
                      'values': {}}
            
            # Store fit results
            phname = str(phmin)+'-'+str(phmax)
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
        
        # Return
        return

    def _save_fits(self):
        """
        Saved phase fit results into a fits file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save phase FITS file')

        # Create FITS file
        self._create_fits()

        # Get file name
        outfile = self['outfile'].filename()

        # Log file name
        self._log_value(gammalib.NORMAL, 'Phase curve file', outfile.url())

        # Save to fits
        self._fits.saveto(outfile, self._clobber())
        
        # Return
        return

    def _save_xml(self):
        """
        Save phase fit results into one XML model each bin
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save phase XML files')

        # Get file name
        outname = self['outfile'].filename()

        # Loop over all phase bins
        for phbin in self._phbins:

            # Set file name for phase bin
            phname        = str(phbin[0])+'-'+str(phbin[1])
            replace_name  = '_phase_'+phname+'.xml'
            outname_phase = (outname.file()).replace('.fits', replace_name)
            outfile       = outname.path() + outname_phase

            # Save XML file
            self._fitmodels[phname].save(outfile)

            # Log file name
            name = 'Phase [%s] file' % phname
            self._log_value(gammalib.NORMAL, name, outfile)
        
        # Return
        return


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
        self._log_observations(gammalib.NORMAL, self._obs, 'Observation')
        
        # Store and clear observations to split them into phase bins
        orig_obs = self._obs.copy()
        
        # Dictionary to save phase fitted models
        self._fitmodels = {}

        # Write header
        self._log_header1(gammalib.TERSE, 'Generate phase curve')

        # Loop over all phases
        for phbin in self._phbins:

            # Write time bin into header
            if self._logTerse():
                self._log.header2('PHASE '+str(phbin[0])+'-'+str(phbin[1]))

            # Select events
            select = ctools.ctselect(orig_obs)
            select['emin'] = self['emin'].real()
            select['emax'] = self['emax'].real()
            select['tmin'] = 'UNDEFINED'
            select['tmax'] = 'UNDEFINED'
            select['rad']  = 'UNDEFINED'
            select['ra']   = 'UNDEFINED'
            select['dec']  = 'UNDEFINED'
            select['expr'] = 'PHASE>'+str(phbin[0])+' && PHASE<'+str(phbin[1])
            select.run()
            
            # Add phase to observation id
            for i in range(0,select.obs().size()):
                oldid = select.obs()[i].id()
                select.obs()[i].id(oldid+'_'+str(phbin[0])+'-'+str(phbin[1]))
            obs = select.obs()
            
            # If a stacked analysis is requested then bin the events
            # and compute the stacked response functions and setup
            # an observation container with a single stacked observation.
            if self._stacked:
                #obs = self._bin_observation(select.obs())
                obs = obsutils.get_stacked_obs(self, select.obs())
    
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
                scaled_norm = model['Normalization'].value()/(phbin[1] - phbin[0])
                model['Normalization'].value(scaled_norm)
            
            # Store fit model
            self._fitmodels[str(phbin[0])+'-'+str(phbin[1])] = like.obs().models().copy()

        # Create FITS file
        self._create_fits()

        # Optionally publish phase curve
        if self['publish'].boolean():
            self.publish()

        # Return
        return

    def save(self):
        """
        Saves results to fits and XML
        """
        # Save results to FITS
        self._save_fits()

        # Save results to XML files (one per phase bin)
        self._save_xml()
        
        # Return
        return

    def publish(self, name=''):
        """
        Publish phase curve

        Parameters
        ----------
        name : str, optional
            Name of phase curve
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Publish phase curve')

        # Continue only if phase curve is valid
        if self._fits.contains('PHASECURVE'):

            # Set default name is user name is empty
            if not name:
                user_name = self._name
            else:
                user_name = name

            # Log file name
            self._log_value(gammalib.NORMAL, 'Phase curve name', user_name)

            # Publish phase curve
            self._fits.publish('PHASECURVE', user_name)

        # Return
        return

    def phasecurve(self):
        """
        Return dictionary with best fit models
        """
        # Return
        return self._fitmodels

        
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csphasecrv(sys.argv)

    # Execute application
    app.execute()
