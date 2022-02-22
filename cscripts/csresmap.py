#! /usr/bin/env python
# ==========================================================================
# Generates a residual map.
#
# Copyright (C) 2014-2022 Michael Mayer
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


# ============== #
# csresmap class #
# ============== #
class csresmap(ctools.csobservation):
    """
    Generates a residual map
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise class members
        self._resmap       = None
        self._use_maps     = False
        self._skip_binning = False

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Initialise some flags
        self._use_maps     = False
        self._skip_binning = False

        # If there are no observations in the observation container then check
        # if the "inobs" parameter is a counts cube. If the "inobs" parameter
        # represents a FITS file and if this FITS file is a counts cube then
        # set self._skip_binning=True
        if self.obs().size() == 0 and self['inobs'].filename() != 'NONE':
            filename = gammalib.GFilename(self['inobs'].filename())
            if filename.is_fits():
                cta = gammalib.GCTAObservation()
                cta.load(filename)
                if cta.eventtype() == 'CountsCube':
                    self._skip_binning = True

        # If we have a counts cube, then ask whether we also have a model cube.
        # If a model cube name is given then set self._use_maps=True
        if self._skip_binning:
            self._use_maps = self['modcube'].is_valid()

        # If not two maps are given, proceed to set up observation
        if not self._use_maps:

            # Set observation if not done before
            if self.obs().size() == 0:
                self._require_inobs('csresmap.get_parameters()')
                self.obs(self._get_observations())

            # If we have exactly one binned CTA observation then signal that
            # the binning can be skipped
            if self.obs().size()         == 1 and \
               self.obs()[0].classname() == 'GCTAObservation' and \
               self.obs()[0].eventtype() == 'CountsCube':
                self._skip_binning = True

            # Set models if we have none
            if self.obs().models().size() == 0:
                self.obs().models(self['inmodel'].filename())

            # If no binning is provided in the observation then query now
            # the counts cube binning parameters
            if not self._skip_binning:
                self['xref'].real()
                self['yref'].real()
                self['coordsys'].string()
                self['proj'].string()
                self['nxpix'].integer()
                self['nypix'].integer()
                self['binsz'].real()
                self['ebinalg'].string()
                if self['ebinalg'].string() == 'FILE':
                    self['ebinfile'].filename()
                else:
                    self['emin'].real()
                    self['emax'].real()
                    self['enumbins'].integer()
                    if self['ebinalg'].string() == 'POW':
                        self['ebingamma'].real()

        # Query parameters
        self['edisp'].boolean()
        self['algorithm'].string()
        self['publish'].boolean()
        self['chatter'].integer()
        self['clobber'].boolean()
        self['debug'].boolean()

        # Query ahead output parameters
        if (self._read_ahead()):
            self['outmap'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return


    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Observation')

        # If a counts and model cube are specified then load them as sky map
        if self._use_maps:
            countmap = gammalib.GSkyMap(self['inobs'].filename())
            modelmap = gammalib.GSkyMap(self['modcube'].filename())

        # ... otherwise build a counts cube and model cube
        else:

            # Do not build counts cube if we have already one in the observation
            # container
            if self._skip_binning:
                cta_counts_cube = gammalib.GCTAEventCube(self.obs()[0].events())

            # ... otherwise generate one now from the event list
            else:

                # Write header
                self._log_header1(gammalib.TERSE, 'Generate binned map (ctbin)')

                # Create counts cube
                cta_counts_cube = obsutils.create_counts_cube(self, self.obs())

            # Assign GCTAEventCube to skymap
            countmap = cta_counts_cube.counts()

            # Write header
            self._log_header1(gammalib.TERSE, 'Generate model map (ctmodel)')

            # Create model map
            model = ctools.ctmodel(self.obs())
            model.cube(cta_counts_cube)
            model['chatter'] = self['chatter'].integer()
            model['clobber'] = self['clobber'].boolean()
            model['debug']   = self['debug'].boolean()
            model['edisp']   = self['edisp'].boolean()
            model.run()

            # Get model map into GSkyMap object
            modelmap = model.cube().counts().copy()

        # Calculate residual maps
        # Note that we need a special
        # construct here to avoid memory leaks. This seems to be a SWIG feature
        # as SWIG creates a new object when calling binning.cube()
        countmap1 = countmap.copy()
        countmap1.stack_maps()
        modelmap.stack_maps()
        self._resmap = obsutils.residuals(self,countmap1,modelmap)

        # Optionally publish map
        if self['publish'].boolean():
            self.publish()

        # Return
        return

    def save(self):
        """
        Save residual map
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save residual map')

        # Get outmap parameter
        outmap = self['outmap'].filename()

        # Continue only filename and residual map are valid
        if self._resmap != None:

            # Log file name
            self._log_value(gammalib.TERSE, 'Residual map file', outmap.url())

            # Save residual map
            self._resmap.save(outmap, self['clobber'].boolean())

            # Stamp residual map
            self._stamp(outmap)

        # Return
        return

    def publish(self, name=''):
        """
        Publish residual map

        Parameters
        ----------
        name : str, optional
            Name of residual map
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Publish residual map')

        # Continue only if residual map is valid
        if self._resmap != None:

            # Set default name is user name is empty
            if not name:
                user_name = self._name()
            else:
                user_name = name

            # Log map name
            self._log_value(gammalib.TERSE, 'Residual map name', user_name)

            # Publish map
            self._resmap.publish(user_name)

        # Return
        return

    def models(self, models):
        """
        Set model
        """
        # Copy models
        self.obs().models(models)

        # Return
        return

    def resmap(self):
        """
        Return residual map

        Returns
        -------
        map : `~gammalib.GSkyMap'
            Residual sky map
        """
        # Return
        return self._resmap


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csresmap(sys.argv)

    # Execute application
    app.execute()
