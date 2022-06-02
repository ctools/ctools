#! /usr/bin/env python
# ==========================================================================
# Puts subset of sources in diffuse model cube
#
# Copyright (C) 2017-2022 Josh Cardenzana
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


# ================= #
# csmodelsois class #
# ================= #
class csmodelsois(ctools.cscript):
    """
    Puts subset of sources in diffuse model cube

    The csmodelsois class puts a subset of sources in a model definition
    XML file into a model cube and generates a new model definition XML file
    in which the subset of sources is replaced by a model of type
    DiffuseMapCube.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor

        Parameters
        ----------
        argv : list of parameters

        Raises
        ------
        TypeError
            An invalid number of command line arguments was provided.
        """
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Initialize parameters
        self._cubegen    = ctools.ctmapcube()
        self._models     = gammalib.GModels()
        self._cubemodels = gammalib.GModels()

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query input parameters
        self['inmodel'].filename()
        self['soilist'].string()

        # Read input parameters
        ptsrcsig = self['ptsrcsig'].real()

        # Get the cube centre coordinates
        self._cubegen['coordsys'].string(self['coordsys'].string())
        if self['coordsys'].string() == 'CEL':
            self._cubegen['xref'].real(self['ra'].real())
            self._cubegen['yref'].real(self['dec'].real())
        else:
            self._cubegen['xref'].real(self['glon'].real())
            self._cubegen['yref'].real(self['glat'].real())

        # Get spatial binning parameters
        self._cubegen['binsz'].real(self['binsz'].real())
        self._cubegen['nxpix'].integer(self['nxpix'].integer())
        self._cubegen['nypix'].integer(self['nypix'].integer())
        self._cubegen['proj'].string(self['proj'].string())

        # Get the energy binning parameters
        self._cubegen['ebinalg'].string(self['ebinalg'].string())
        if self['ebinalg'].string() == 'FILE':
            self._cubegen['ebinfile'].filename(self['ebinfile'].filename())
        else:
            self._cubegen['emin'].real(self['emin'].real())
            self._cubegen['emax'].real(self['emax'].real())
            self._cubegen['enumbins'].integer(self['enumbins'].integer())

        # Set point source significance
        self._cubegen['ptsrcsig'].real(ptsrcsig)

        # Read optionally output cube filenames
        if self._read_ahead():
            self._cubegen['outcube'].filename(self['outcube'].filename())
            self['outmodel'].query()

        # Write parameters into logger
        self._log_parameters(gammalib.TERSE)

        return

    def _gen_cubemodel(self):
        """
        Generates a binned model from the model file using ctmapcube
        """
        # Get a list of sources we DONT want to put into the cube
        sources = self['soilist'].string().split(',')

        # Load the model into a model container
        if (self._models.size() == 0):
            self._models = gammalib.GModels(self['inmodel'].filename())

        # Store a copy of the models
        self._cubemodels = self._models.copy()

        # Loop through models and pull out all models not in the list
        for model in self._models:

            # If model is not a sky model then continue
            if model.classname() != 'GModelSky':
                self._cubemodels.remove(model.name())

            # ... otherwise, if model is in list of sources the remove it
            # from the container
            elif model.name() in sources:
                self._cubemodels.remove(model.name())

        # Log the number of models to be put into the cube
        self._log_value(self['chatter'].integer(), 'Numner of cube models',
                        self._cubemodels.size())

        # Return
        return


    # Public methods
    def process(self):
        """
        Implements the actual bulk of the script's tasks
        """
        # Get parameters
        self._get_parameters()

        # Set the models to put into the cube
        self._gen_cubemodel()

        # Only run if there still models to be filled into the cube
        if not self._cubemodels.is_empty():

            # Log action
            self._log_string(gammalib.NORMAL, 'Generating model cube')

            # Fill the cube models into the ctmapcube object
            self._cubegen.models(self._cubemodels)

            # Run the map generating method
            self._cubegen.run()

        else:
            self._log_string(gammalib.NORMAL,
                             'List of models is empty, nothing will be done')

        # Return
        return

    def save(self):
        """
        Save the cube and optionally the updated XML model file
        """
        # Save the generated cube if the cube and filename are not empty
        if ((not self._cubegen.mapcube().cube().is_empty()) and
            self['outcube'].is_valid()):

            # Save cube
            self._cubegen.mapcube().save(self['outcube'].filename(),
                                         self['clobber'].boolean())

            # Stamp cube
            self._stamp(self['outcube'].filename())

        # If requested, save the updated list of models
        if self['outmodel'].is_valid():

            # Generate a list of models that will be output to file
            outmodels = gammalib.GModels(self._models)

            # Remove all models used in the generated cube
            for model in self._cubemodels:
                outmodels.remove(model.name())

            # Generate the actual cube model
            spat     = gammalib.GModelSpatialDiffuseCube(self['outcube'].filename())
            spec     = gammalib.GModelSpectralConst()
            newmodel = gammalib.GModelSky(spat, spec)

            # Give the model a default name
            newmodel.name(self.cubemodelname())

            # Now append the model to the list of output models
            outmodels.append(newmodel)

            # Save the model list
            outmodels.save(self['outmodel'].filename())

        # Return
        return

    def mapcube(self):
        """
        Return the mapcube generated by the underlying 'ctmapcube' object
        """
        # Return
        return self._cubegen.mapcube()

    def cubemodelname(self):
        """
        Return the name to be given to the cube model in the output XML file
        """
        # Return
        return 'FixedSourcesCube'

    def models(self, models):
        """
        Set the models

        Parameters
        -------
        models : `~gammalib.GModels`
            Model container
        """
        # Set models
        self._models = models

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
