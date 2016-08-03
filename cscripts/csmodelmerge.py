#!/usr/bin/env python
# ==========================================================================
# Merge model definition XML files
#
# Copyright (C) 2015-2016 Michael Mayer
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
import os
import sys
import glob
import gammalib
import ctools


# ================== #
# csmodelmerge class #
# ================== #
class csmodelmerge(ctools.cscript):
    """
    Merge model definition XML files

    An arbitrary number of model definition XML files will be merged into
    a single model definition XML file.    
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = 'csmodelmerge'
        self._version = '1.2.0'

        # Initialise class members
        self._files      = None
        self._models     = gammalib.GModels()

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Get input models string
        inmodels = self['inmodels'].string()
        
        # Handle ASCII files. If the file names given in the ASCII are
        # relative filenames it is assumed that the filename is given
        # relative to the location of the file.
        if '@' == inmodels[0]:
            filename    = inmodels.replace('@','')
            self._files = open(filename).read().splitlines()
            dirname     = os.path.dirname(filename)
            files       = []
            for f in self._files:
                if f[0] != '/':
                    fname = dirname + '/' + f
                else:
                    fname = f
                files.append(fname)
            self._files = files
    
        # Handle wild card strings
        elif '*' in inmodels:
            self._files = glob.glob(inmodels)
        
        # Handle space separated list
        elif ' ' in inmodels:
            self._files = inmodels.split(' ')
        
        # Handle semi-colon separated list
        elif ';' in inmodels:
            self._files = inmodels.split(';')
            
        # Throw exception if input models cannot be decoded
        else:
            msg = 'Parameter "inmodels" must contain either an @ASCII '\
                  'file, a semi-colon-separated or whitespace-separated '\
                  'list of files or a wildcard string.'
            raise RuntimeError(msg)
        
        # Read ahead output filename
        if self._read_ahead():
            self['outmodel'].filename()
            
        # Get clobber parameter
        self._clobber = self['clobber'].boolean()

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
        
        # Write header
        if self._logTerse():
            self._log('\n')
            self._log.header1('Merge models')

        # Initialise model container
        self._models = gammalib.GModels()
        
        # Loop over model files
        for f in self._files:

            # Construct container from XML file
            models = gammalib.GModels(f)
            
            # Log number of models to add
            nmodels = models.size()
            if nmodels == 0:
                name = 'Add no model from file'
            elif nmodels == 1:
                name = 'Add 1 model from file'
            else:
                name = 'Add %d models from file' % nmodels
            self._log_value(gammalib.TERSE, name, f)
            
            # Extend model container by adding all models in the model file
            self._models.extend(models)
                
        # Log total number of models
        self._log_value(gammalib.TERSE, 'Models after merging',
                                        self._models.size())

        # Return
        return

    def save(self):
        """ 
        Save model definition XML file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save models')

        # Get output filename in case it was not read ahead
        outmodel = self['outmodel'].filename()
        
        # If file exists and clobber flag is false then raise an exception
        if outmodel.exists() and not self._clobber:
            msg = ('Cannot save "'+outmodel.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)
        
        else:
            
            # Log filename
            self._log_value(gammalib.NORMAL, 'Model definition XML file',
                                             outmodel.url())
    
            # Save models
            self._models.save(outmodel)
        
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

        # Save ds9 file if required
        self.save()

        # Return
        return    
        

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csmodelmerge(sys.argv)
    
    # Execute application
    app.execute()
