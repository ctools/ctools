#! /usr/bin/env python
# ==========================================================================
# Convolve models with response for COMPTEL observations
#
# Copyright (C) 2022 Juergen Knoedlseder
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
import glob
import os
import gammalib
import ctools


# ================ #
# comobsconv class #
# ================ #
class comobsconv(ctools.csobservation):
    """
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Set observation if not done before
        if self.obs().is_empty():
            self.obs().load(self['inobs'].filename())

        # Set models if we have none
        if self.obs().models().is_empty():
            self.obs().models(self['inmodel'].filename())

        # Query parameters
        self['outfolder'].string()
        self['filetype'].string()

        # Query ahead output model filename
        if self._read_ahead():
            self['outobs'].filename()

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

        # Set parameters
        outfolder = self['outfolder'].string()
        filetype  = self['filetype'].string()

        # Create cache file directory
        try:
            os.makedirs(gammalib.expand_env(outfolder))
        except OSError:
            pass

        # Write input observation container into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Input observation')

        # Write input models into logger
        self._log_models(gammalib.VERBOSE, self.obs().models(), 'Input model')

        # Log header
        self._log_header1(gammalib.NORMAL, 'Set response cache names')

        # Loop over all observations
        for obs in self.obs():

            # Write header
            self._log_header3(gammalib.NORMAL, self._get_obs_header(obs))

            # Get DRB name
            drbname = obs.drbname().url()

            # Derive response cache name from DRB name
            rspname = '%s/%s' % (outfolder, os.path.basename(drbname.replace('drb', filetype)))

            # Set response cache name
            obs.rspname(rspname)

            # Log DRE and cache name
            self._log_value(gammalib.NORMAL, 'DRB name', drbname)
            self._log_value(gammalib.NORMAL, 'Response cache name', rspname)

        # Log header
        self._log_header1(gammalib.NORMAL, 'Convolve models with response')

        # Compute response by calling evaluation method
        self.obs().eval()

        # Return
        return

    def save(self):
        """ 
        Save observation definition file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save observations')

        # Get output filename
        outobs = self['outobs'].filename()

        # If file exists and clobber flag is false then raise an exception
        if outobs.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outobs.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:

            # Save observation definition file
            self.obs().save(outobs)
            self._log_value(gammalib.NORMAL, 'Obs. definition XML file', outobs.url())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comobsconv(sys.argv)

    # Execute application
    app.execute()
