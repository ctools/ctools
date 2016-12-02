#!/usr/bin/env python
# ==========================================================================
# Detects sources in a sky map
#
# Copyright (C) 2016 Juergen Knoedlseder
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
import math
import gammalib
import ctools
from cscripts import modutils


# ================= #
# cssrcdetect class #
# ================= #
class cssrcdetect(ctools.cscript):
    """
    Detects sources in a sky map
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'cssrcdetect'
        self._version = '1.2.0'

        # Set protected members
        self._models = gammalib.GModels()

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query input parameters
        self['inmap'].filename()

        # Query further parameters
        self['srcmodel'].string()
        self['bkgmodel'].string()
        self['threshold'].real()
        self['maxsrcs'].integer()
        self['exclrad'].real()
        self['fit_pos'].boolean()
        self['fit_shape'].boolean()
        
        # Query ahead output model filename
        if self._read_ahead():
            self['outmodel'].filename()
            self['outds9file'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _detect_sources(self, counts):
        """
        Detect sources in counts map

        Parameters
        ----------
        counts : `~gammalib.GSkyMap()`
            Counts map
        """
        # Source detection iterations
        for i in range(self['maxsrcs'].integer()):

            # Write header
            self._log_header3(gammalib.NORMAL, 'Iteration '+str(i+1))

            # Get map moments
            mean, std = self._map_moments(counts)

            # Log map moments
            self._log_value(gammalib.NORMAL, 'Map mean', mean)
            self._log_value(gammalib.NORMAL, 'Map standard deviation', std)
 
            # Compute threshold
            threshold = mean + self['threshold'].real() * std

            # Log map threshold
            self._log_value(gammalib.NORMAL, 'Map threshold', threshold)

            # Get maximum value and corresponding sky direction
            value, dir = self._find_maximum(counts, threshold=threshold)

            # If maximum found then log maximum and add model
            if dir is not None:

                # Set source name
                name = 'Src%3.3d' % (i+1)

                # Log maximum
                self._log_value(gammalib.NORMAL, 'Map maximum', str(value))
                self._log_value(gammalib.NORMAL, name+' position', str(dir))
            
                # Add model
                self._add_model(dir, name)

                # Remove maximum from map
                counts = self._remove_maximum(counts, dir, mean,
                                              radius=self['exclrad'].real())

            # ... otherwise log that no maximum was found and break iterations
            else:

                # Log than no maximum was found
                self._log_value(gammalib.NORMAL, 'Map maximum',
                                'None above threshold')

                # Break
                break

        # Return
        return

    def _find_maximum(self, map, threshold=0.0):
        """
        Find maximum in a sky map

        Parameters
        ----------
        map : `~gammalib.GSkyMap()`
            Sky map
        threshold : float, optional
            Threshold for maximum value

        Returns
        -------
        value, dir : tuple of float and `~gammalib.GSkyDir()`
            Maximum sky map value and corresponding sky direction
        """
        # Initialise maximum pixel value and sky direction
        value = threshold
        dir   = None

        # Loop over all pixels and find maximum
        for i in range(map.npix()):
            if map[i] > value:
                value = map[i]
                dir   = map.inx2dir(i)

        # Return sky direction of maximum
        return value, dir

    def _remove_maximum(self, map, dir, value=0.0, radius=0.1):
        """
        Remove maximum from sky map by replacing pixels with a given value

        Parameters
        ----------
        map : `~gammalib.GSkyMap()`
            Sky map
        dir : `~gammalib.GSkyDir()`
            Sky direction of maximum
        value : float, optional
            Replacement value
        radius : float, optional
            Radius within which pixel values are replaced

        Returns
        -------
        map : `~gammalib.GSkyMap()`
            Sky map with maximum removed
        """
        # Copy skymap
        map_copy = map.copy()
        
        # Loop over all pixels and find maximum
        for i in range(map_copy.npix()):
            mapdir = map_copy.inx2dir(i)
            if mapdir.dist_deg(dir) < radius:
                map_copy[i] = value

        # Return copy of map
        return map_copy

    def _map_moments(self, map):
        """
        Determine moments of sky map pixels

        Parameters
        ----------
        map : `~gammalib.GSkyMap()`
            Sky map

        Returns
        -------
        mean, std : tuple of float
            Mean and standard deviation of pixel values
        """
        # Initialise mean and standard deviation
        mean = 0.0
        std  = 0.0

        # Loop over all pixels
        for i in range(map.npix()):
            mean += map[i]
            std  += map[i]*map[i]

        # Compute mean and standard deviation
        mean /= float(map.npix())
        std  /= float(map.npix())
        std  -= mean
        std   = math.sqrt(std)

        # Return mean and standard deviation
        return mean, std

    def _add_model(self, dir, name):
        """
        Add model to model container

        Parameters
        ----------
        dir : `~gammalib.GSkyDir()`
            Sky direction of model
        name : str
            Model name
        """
        # Set point source model
        model = self._set_ptsrc(dir)

        # Set model name
        model.name(name)

        # Append model to container
        self._models.append(model)

        # Return
        return

    def _set_bkg(self, type):
        """
        Set background model

        Parameters
        ----------
        type : str
            Model type ('IRF', 'AEFF' or 'CUBE')

        Returns
        -------
        model : `~gammalib.GModelData()`
            Background model
        """
        # Set spectral component
        spectral = gammalib.GModelSpectralPlaw(1.0, 0.0,
                                               gammalib.GEnergy(1.0, 'TeV'))

        # Set background model
        if type == 'IRF':
            model = gammalib.GCTAModelIrfBackground(spectral)
        elif type == 'AEFF':
            model = gammalib.GCTAModelAeffBackground(spectral)
        elif type == 'CUBE':
            model = gammalib.GCTAModelCubeBackground(spectral)
        else:
            model = None

        # Set background name
        if model is not None:
            model.name('Background')

        # Return model
        return model

    def _set_ptsrc(self, dir):
        """
        Set point source model

        Parameters
        ----------
        dir : `~gammalib.GSkyDir()`
            Sky direction of model

        Returns
        -------
        model : `~gammalib.GModelSky()`
            Point source model
        """
        # Set spatial component
        spatial = gammalib.GModelSpatialPointSource(dir)

        # Get fit parameters
        fit_pos   = self['fit_pos'].boolean()
        fit_shape = self['fit_shape'].boolean()

        # Loop over all parameters
        for par in spatial:

            # Handle position parameters
            if par.name() == 'RA' or par.name() == 'DEC':
                if fit_pos:
                    par.free()
                else:
                    par.fix()

        # Set spectral component (powerlaw with 1% of Crab)
        spectral = gammalib.GModelSpectralPlaw(5.7e-18, -2.48,
                                               gammalib.GEnergy(0.3, 'TeV'))

        # Set sky model
        model = gammalib.GModelSky(spatial, spectral)

        # Return model
        return model


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Clear model container
        self._models.clear()

        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Write header
        self._log_header1(gammalib.NORMAL, 'Source detection')

        # Get skymap filename
        inmap = self['inmap'].filename()

        # Open sky map file
        fits = gammalib.GFits(inmap)

        # Extract primary extension as counts map
        counts = gammalib.GSkyMap(fits.image(0))

        # Detect sources
        self._detect_sources(counts)

        # Close sky map file
        fits.close()

        # Write detected sources as models into logger
        self._log_models(gammalib.NORMAL, self._models, 'Detected source model')

        # Write header
        self._log_header1(gammalib.NORMAL, 'Add background model')

        # Add background model if requested
        bkgmodel = self._set_bkg(self['bkgmodel'].string())
        if bkgmodel is not None:

            # Append model
            self._models.append(bkgmodel)

            # Log model
            self._log_string(gammalib.NORMAL, str(bkgmodel))

        # ... otherwise notify that no background model is added
        else:
            self._log_value(gammalib.NORMAL, 'Background model', 'None')

        # Return
        return

    def save(self):
        """ 
        Save sources into model definition XML file and ds9 file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save sources')

        # Get output filenames
        outmodel   = self['outmodel'].filename()
        outds9file = self['outds9file'].filename()

        # If model definition XML file exists and clobber flag is false then
        # raise an exception
        if outmodel.exists() and not self._clobber:
            msg = ('Cannot save "'+outmodel.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # If ds9 file exists and clobber flag is false then raise an exception
        elif outds9file.exists() and not self._clobber:
            msg = ('Cannot save "'+outds9file.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:
            # Log model definition XML filename
            self._log_value(gammalib.NORMAL, 'Model definition XML file',
                                             outmodel.url())

            # Save models
            self._models.save(outmodel)

            # Log ds9 filename
            self._log_value(gammalib.NORMAL, 'DS9 region file',
                                             outds9file.url())

            # Save models as ds9 region file
            modutils.models2ds9file(self._models, outds9file.url())

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

        # Save observation definition file
        self.save()

        # Return
        return    


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cssrcdetect(sys.argv)

    # Execute application
    app.execute()
