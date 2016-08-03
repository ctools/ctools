#!/usr/bin/env python
# ==========================================================================
# Find observations from an IACT data store
#
# Copyright (C) 2016 Michael Mayer
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
import os
import json
import gammalib
import ctools


# =============== #
# csfindobs class #
# =============== #
class csfindobs(ctools.cscript):
    """
    Find observations from an IACT data store
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name         = 'csfindobs'
        self._version      = '1.2.0'
        self._datapath     = os.getenv('VHEFITS','')
        self._prodname     = ''
        self._select_radec = True
        self._radius       = 0.0
        self._ra           = 0.0
        self._dec          = 0.0
        self._obs_index    = ''
        self._runs         = []
        
        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return
    

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters
        if self._datapath == '':
            self._datapath = self['datapath'].string()
        
        # Expand environment
        self._datapath = gammalib.expand_env(self._datapath)
        
        # Get production name
        self._prodname = self['prodname'].string()
        
        # Master index file name
        master_indx = self['master_indx'].string()
        
        # Initialise flag if spatial selection is required
        self._select_radec = True
        
        # Initialise invalid radius
        self._radius = 0.0
        
        # Check for validity of spatial parameters
        if (self['ra'].is_valid() and
            self['dec'].is_valid() and 
            self['rad'].is_valid()):
            
            # Read spatial parameters
            self._ra     = self['ra'].real()
            self._dec    = self['dec'].real()
            self._radius = self['rad'].real()

        # ... otherwise signal that there are no spatial parameters for
        # selection
        else:
            self._select_radec = False
        
        # Check Radius for validity
        if self._radius <= 0.0:
            self._select_radec = False
            
        # Query other parameters
        self['min_qual'].integer()
        self['expression'].string()

        # Read ahead output parameters
        if self._read_ahead():
            self['outfile'].filename()
        
        # Set filename of JSON master file and raise an exception if the file
        # does not exist
        master_file = os.path.join(self._datapath, master_indx)
        if not os.path.isfile(master_file):
            msg = ('FITS data store not available. No master index file found '
                   'at "%s". Make sure the file is copied from the server and '
                   'your datapath is set correctly.' % master_file)
            raise RuntimeError(msg)

        # Open and load JSON master file. If the "dataset" key is not available
        # then raise an exception
        json_data = open(master_file).read()
        data      = json.loads(json_data)
        if not 'datasets' in data:
            msg = ('Key "datasets" not available in master index file.')
            raise RuntimeError(msg)

        # Get configurations from JSON master file
        configs = data['datasets']

        # Initialise obs index file
        self._obs_index = ''

        # Get name of observation index file
        for config in configs:
            if self._prodname == config['name']:
                self._obs_index = str(os.path.join(self._datapath,
                                                   config['obsindx']))
                break

        # If the observation index file name is empty then raise an exception
        if self._obs_index == '':
            msg = ('FITS data store "%s" not available. Run csiactdata to get '
                   'a list of available storage names.' % self._prodname)
            raise RuntimeError(msg)

        # If the observation index file is not a FITS file then raise an
        # exception
        filename = gammalib.GFilename(self._obs_index+'[OBS_INDEX]')
        if not filename.is_fits():
            msg = ('Observation index file "%s[OBS_INDEX]" for FITS data store '
                   '"%s" not available. Check your master index file or run '
                   'csiactdata to get a list of available storage names.' %
                   (self._obs_index, self._prodname))
            raise RuntimeError(msg)

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

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

        # Write header into logger
        self._log_header1(gammalib.TERSE, 'Find observations')
        
        # Initialise run list
        self._runs = []
        
        # Initialise selection expression
        expr = ''
        
        # If spatial selection is requested then add an angular separation
        # expression to the expression string
        if self._select_radec:
            expr += 'ANGSEP('+str(self._ra)+','+str(self._dec)+ \
                    ',RA_PNT,DEC_PNT)<='+str(self._radius)
        
        # Add '&&' connector if expression is not empty
        if len(expr):
            expr += '&&'
        
        # Add always quality expression
        expr += 'QUALITY<='+str(self['min_qual'].integer())
        
        # Add user expression is one has been specified
        expression = self['expression'].string()    
        if (expression != 'NONE' and expression != 'INDEF' and
            len(expression) > 0):
            expr += '&&' + expression   

        # Write the expression into the logger
        self._log_value(gammalib.NORMAL, 'Expression', expr)

        # Set filename including the selection expression
        filename = self._obs_index+'[OBS_INDEX]['+expr+']'
        
        # Open observation index FITS file
        fits = gammalib.GFits(filename)

        # If there are entries in the observation index table then append the
        # observations identifiers to the run list
        obs_index = fits['OBS_INDEX']
        if obs_index.nrows() > 0:
            for i in range(obs_index.nrows()):
                self._runs.append(obs_index['OBS_ID'][i])

        # Write the number of observations into the logger
        self._log_value(gammalib.TERSE, 'Observations', len(self._runs))

        # Write the observation identifiers into the logger
        if len(self._runs) > 0:
            for i, run in enumerate(self._runs):
                self._log_value(gammalib.NORMAL, 'Observation %d' % (i+1), run)

        # Close FITS file
        fits.close()

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
        
        # Save residual map
        self.save()
        
        # Return
        return

    def save(self):
        """
        Save runlist
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save runlist')

        # Get filename
        outfile = self['outfile'].filename()

        # If the runlist file exists but the clobber flag is not set then signal
        # that runlist file exists already
        if outfile.exists() and not self._clobber():
            msg = ('File "'+outfile+'" exists already, runlist not saved. '
                   'Set clobber=yes to overwrite the file.\n')
            self._log_string(gammalib.TERSE, msg)

        # ... otherwise save the runlist into a file
        else:

            # Log file name
            self._log_value(gammalib.NORMAL, 'Runlist file', outfile.url())

            # Write all runs to file
            f = open(outfile.url(),'w')
            for run in self._runs:
                f.write(str(run)+' \n')
            f.close()

        # Return
        return

    def runs(self):
        """ 
        Return run list

        Returns
        -------
        runs : list of str
            A list of runs
        """
        # Return
        return self._runs
    

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csfindobs(sys.argv)
    
    # Execute application
    app.execute()
