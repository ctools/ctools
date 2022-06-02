#! /usr/bin/env python
# ==========================================================================
# Simulate COMPTEL observations
#
# Copyright (C) 2021-2022 Juergen Knoedlseder
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
import math
import gammalib
import ctools


# =============== #
# comobssim class #
# =============== #
class comobssim(ctools.csobservation):
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
        self['suffix'].string()
        self['seed'].integer()
        self['add'].boolean()

        # Query ahead output skymap filename
        if self._read_ahead():
            self['outfolder'].string()
            self['outobs'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _simulate_events(self, dre, drm, ran):
        """
        Simulate events

        Parameters
        ----------
        dre : `~gammalib.GCOMDri`
            Event cube
        drm : `~gammalib.GCOMDri`
            Model cube
        ran : `~gammalib.GRan`
            Random number generator

        Returns
        -------
        dre : `~gammalib.GCOMDri`
            Event cube including simulated events
        """
        # Get simulation mode
        add = self['add'].boolean()

        # Loop over all dataspace pixels
        for i in range(drm.size()):

            # Simulate events
            events = ran.poisson(drm[i])

            # Set or add events
            if add:
                dre[i] += events
            else:
                dre[i] = events

        # Return DRE
        return dre


    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Log header
        self._log_header1(gammalib.NORMAL, 'Input observations')

        # Log input observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

        # Write header
        self._log_header1(gammalib.NORMAL, 'Input models')

        # Log models
        self._log_string(gammalib.NORMAL, str(self.obs().models()))

        # Write header
        self._log_header1(gammalib.NORMAL, 'Simulate events')

        # Set random number generator seed
        seed = self['seed'].integer()
        ran  = gammalib.GRan(seed)

        # Get DRE file suffix
        suffix = self['suffix'].string()

        # Write seed
        self._log_value(gammalib.NORMAL, 'Seed', seed)

        # Loop over all input observations
        for obs in self.obs():

            # Write header
            self._log_header3(gammalib.NORMAL, self._get_obs_header(obs))

            # Get DRE
            dre = obs.events().dre()

            # Compute DRM
            drm = obs.drm(self.obs().models())

            # Simulate events
            dre = self._simulate_events(dre, drm, ran)

            # Set DRE
            events = gammalib.GCOMEventCube(dre)
            obs.events(events)

            # Set DRE filename
            dreid   = 'dre_sim-seed%6.6d' % (seed)
            if suffix != '':
                dreid += '%s_' % (suffix)
            drename = os.path.basename(obs.drename().url())
            drename = '%s/%s' % (self['outfolder'].string(),
                                 drename.replace('dre', dreid))
            obs.drename(drename)

        # Write header
        self._log_header1(gammalib.NORMAL, 'Simulated observations')

        # Log simulated observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

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

            # Get outfolder
            outfolder = self['outfolder'].string()

            # Create outfolder directory
            try:
                os.makedirs(gammalib.expand_env(outfolder))
            except OSError:
                pass

            # Log filename
            self._log_value(gammalib.NORMAL, 'Obs. definition XML file',
                                             outobs.url())

            # Save observation definition XML file
            self.obs().save(outobs)

            # Save simulated DREs
            for obs in self.obs():

                # Log filename
                self._log_value(gammalib.NORMAL, 'DRE file',
                                                 obs.drename().url())

                # Save DRE file
                dre = obs.events().dre()
                dre.save(obs.drename().url(), self['clobber'].boolean())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comobssim(sys.argv)

    # Execute application
    app.execute()
