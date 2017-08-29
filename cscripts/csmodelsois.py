#! /usr/bin/env python
# ==========================================================================
# Allows generation of DiffuseMapCube model from subset of sources in a
# model definition XML file
#
# Copyright (C) 2017 Josh Cardenzana
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

#=============#
# csmodelsois #
#=============#
class csmodelsois(ctools.cscript):
    """
    Generates a DiffuseMapCube model from a subset of sources in a model
    definition XML file. All background models will be ignored!
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        
        Parameters
        ----------
        argv : list of parameters
        
        Raises
        ------
        TypeError
        An invalid number of command line arguments was provided.
        """
        # Set name and version
        self._name    = 'csmodelsois'
        self._version = ctools.__version__
        
        # Initialize parameters
        self._init_members()

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)
        
        # Return
        return


    def _init_members(self):
        """
        Initialize the parameters for this object
        """
        self.cubegen    = ctools.ctmapcube()    # Implements ctmapcube for generate the cube
        self.m_models   = gammalib.GModels()    # Stores the full list of input models
        self.cubemodels = gammalib.GModels()    # Stores the list of models to be written to file
        return


    def _set_cube_parameters(self):
        """
        Fill in the parameters associated with this object
        """
        
        # Get the cube centre coordinates
        self.cubegen['coordsys'].string( self['coordsys'].string() )
        if self['coordsys'].string() == 'CEL':
            # Center in celestial coordinates
            self.cubegen['xref'].real( self['ra'].real() )
            self.cubegen['yref'].real( self['dec'].real() )
        else:
            # Otherwise assume center in galactic coordinates
            self.cubegen['xref'].real( self['glon'].real() )
            self.cubegen['yref'].real( self['glat'].real() )

        # Get spatial binning parameters
        self.cubegen['binsz'].real( self['binsz'].real() )
        self.cubegen['nxpix'].integer( self['nxpix'].integer() )
        self.cubegen['nypix'].integer( self['nypix'].integer() )
        self.cubegen['proj'].string( self['proj'].string() )
        
        # Get the energy binning parameters
        self.cubegen['ebinalg'].string( self['ebinalg'].string() )
        if self['ebinalg'].string() == "FILE":
            # Read energy binning from a file
            self.cubegen['ebinfile'].filename( self['ebinfile'].filename() )
        else:
            # Create the energy binning information
            self.cubegen['emin'].real( self['emin'].real() )
            self.cubegen['emax'].real( self['emax'].real() )
            self.cubegen['enumbins'].integer( self['enumbins'].integer() )

        # Get remaining hidden parameters
        self.cubegen['ptsrcsig'].real( self['ptsrcsig'].real() )
        
        # Read optionally output cube filenames
        if self._read_ahead():
            self.cubegen['outcube'].filename( self['outcube'].filename() )
            self['soilist'].string()
            self['outmodel'].filename()

        # Write parameters into logger
        self._log_parameters(gammalib.TERSE);

        return


    def _gen_cubemodel(self):
        """
        Generates a binned model from the model file using ctmapcube
        returns: updated model filename
        """
        # Get a list of sources we DONT want to put into the cube
        sources = self['soilist'].string().split(',')

        # Load the model into a model container
        if (self.m_models.size() == 0):
            self.m_models = gammalib.GModels(self['inmodel'].filename())
        
        self.cubemodels = self.m_models.copy()

        # Loop through models and pull out all models not in the list
        models2remove = []
        for model in self.m_models:
            # Convert the model to a GModelSky
            mod = model.classname()

            # Check the model
            if mod != "GModelSky":
                # Model is most likely a background model, so continue
                self.cubemodels.remove(model.name())
            elif model.name() in sources:
                # Model is a source of interest, so remove it from list of models
                # to be used in generating the cube
                self.cubemodels.remove(model.name())

        # Log the number of models to be put into the cube
        self._log_value(self['chatter'].integer(), 'Cube models', self.cubemodels.size())

        return


    def mapcube(self):
        """
        Return the mapcube generated by the underlying 'ctmapcube' object
        """
        return self.cubegen.mapcube()


    def cubemodelname(self):
        """
        Return the name to be given to the cube model in the output XML file
        """
        return "FixedSourcesCube"


    def models(self, models:gammalib.GModels):
        """
        Set the models for analysis
        """
        self.m_models = models
        return


    def run(self):
        """
        Implements the actual bulk of the script's tasks
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._set_cube_parameters()

        # Set the models to put into the cube
        self._gen_cubemodel()

        # Only run if there still models to be filled into the cube
        if not self.cubemodels.is_empty():
            self._log_string(gammalib.NORMAL, 'Generting model cube')

            # Fill the cube models into the ctmapcube object
            self.cubegen.models( self.cubemodels )

            # Run the map generating method
            self.cubegen.run()

        else:
            self._log_string(gammalib.NORMAL, 'List of models is empty, nothing will be done')

        return


    def save(self):
        """
        Save the cube and optionally the updated XML model file
        """
        #cube = gammalib.GModelSpatialDiffuseCube(self.cubegen.mapcube())

        # Save the generated cube if the cube and filename are not empty
        if ((not self.cubegen.mapcube().cube().is_empty()) and
            self._is_valid_filename(self['outcube'].filename())):

            self.cubegen.mapcube().save( self['outcube'].filename(), self['clobber'].boolean() )
        
        # If requested, save the updated list of models
        if self._is_valid_filename( self['outmodel'].filename() ):
            # Generate a list of models that will be output to file
            outmodels = gammalib.GModels(self.m_models)
            
            # Remove all models used in the generated cube
            for model in self.cubemodels:
                outmodels.remove(model.name())

            # Generate the actual cube model
            spat     = gammalib.GModelSpatialDiffuseCube( self['outcube'].filename() )
            spec     = gammalib.GModelSpectralConst()
            newmodel = gammalib.GModelSky(spat, spec)

            # Give the model a default name
            newmodel.name( self.cubemodelname() )

            # Now append the model to the list of output models
            outmodels.append(newmodel)

            # Save the model list
            outmodels.save( self['outmodel'].filename() )

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
        self.save()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    
    # Create instance of application
    app = csmodelsois(sys.argv)
    
    # Execute application
    app.execute()
