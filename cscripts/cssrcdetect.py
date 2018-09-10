#!/usr/bin/env python
# ==========================================================================
# Detects sources in a sky map
#
# Copyright (C) 2016-2018 Juergen Knoedlseder
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
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Set protected members
        self._models = gammalib.GModels()

        # Initialise sky map from constructor arguments
        if len(argv) > 0 and isinstance(argv[0], gammalib.GSkyMap):
            self._map = argv[0]
            argv      = argv[1:]
        else:
            self._map = gammalib.GSkyMap()

        # Initialise other members
        self._map_dirs = []

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query input parameters if sky map is empty
        if self._map.is_empty():
            self['inmap'].filename()

        # Query further parameters
        self['srcmodel'].string()
        self['bkgmodel'].string()
        self['threshold'].real()
        self['maxsrcs'].integer()
        self['exclrad'].real()
        self['fit_pos'].boolean()
        self['fit_shape'].boolean()

        # Query the smoothing parameters
        self['smoothkernel'].string()
        if self['smoothkernel'].string().upper() != 'NONE':
            self['smoothparam'].real()

        # Query ahead output model filename
        if self._read_ahead():
            self['outmodel'].filename()
            self['outds9file'].filename()
            self['momentradius'].real()

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
            mean, std = self._map_moments(counts, self['momentradius'].real())
 
            # Compute threshold
            sigmap = (counts - mean)/std
            
            # Get maximum value and corresponding sky direction
            value, pos = self._find_maximum(sigmap)

            # If maximum found then log maximum and add model
            if pos is not None:

                # Set source name
                name = 'Src%3.3d' % (i+1)

                # Log maximum
                self._log_value(gammalib.NORMAL, 'Map maximum', str(value))
                self._log_value(gammalib.NORMAL, name+' position', str(pos))
            
                # Add model
                self._add_model(pos, name)

                # Remove maximum from map
                counts = self._remove_maximum(counts, pos, mean(pos),
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

    def _find_maximum(self, skymap):
        """
        Find maximum in a sky map

        Parameters
        ----------
        skymap : `~gammalib.GSkyMap()`
            Sky map

        Returns
        -------
        value, pos : tuple of float and `~gammalib.GSkyDir()`
            Maximum sky map value and corresponding sky direction
        """
        # Initialise maximum pixel value and sky direction
        value = self['threshold'].real()
        pos   = None

        # Loop over all pixels and find maximum
        for i in range(skymap.npix()):
            if skymap[i] > value:
                value = skymap[i]
                pos   = self._map_dirs[i]

        # Return sky direction of maximum
        return value, pos

    def _remove_maximum(self, skymap, pos, value=0.0, radius=0.1):
        """
        Remove maximum from sky map by replacing pixels with a given value

        Parameters
        ----------
        skymap : `~gammalib.GSkyMap()`
            Sky map
        pos : `~gammalib.GSkyDir()`
            Sky direction of maximum
        value : float, optional
            Replacement value
        radius : float, optional
            Radius within which pixel values are replaced

        Returns
        -------
        skymap : `~gammalib.GSkyMap()`
            Sky map with maximum removed
        """
        # Copy skymap
        skymap_copy = skymap.copy()
        
        # Cache the cosine of the radius
        cos_radius = math.cos(math.radians(radius))

        # Loop over all pixels and find maximum
        for i in range(skymap_copy.npix()):
            skymap_dir = self._map_dirs[i]
            if skymap_dir.cos_dist(pos) > cos_radius:
                skymap_copy[i] = value

        # Return copy of map
        return skymap_copy

    def _map_moments(self, skymap, radius):
        """
        Determine moments of sky map pixels

        Parameters
        ----------
        skymap : `~gammalib.GSkyMap()`
            Sky map
        radius : float
            radius (deg) for pixel consideration

        Returns
        -------
        mean, std : tuple of GSkyMap
            Mean and standard deviation of pixel values within a given radius
        """
        # Copy the input skymap
        mean = gammalib.GSkyMap(skymap)
        std  = gammalib.GSkyMap(skymap)
        std *= std

        # Convolve by disk to get bin-by-bin mean
        mean.smooth('DISK', radius)
        std.smooth('DISK', radius)

        # Compute the standard deviation for each pixel
        std = std - (mean*mean)
        std = std.sqrt()

        # Return mean and standard deviation
        return mean, std

    def _add_model(self, pos, name):
        """
        Add model to model container

        Parameters
        ----------
        pos : `~gammalib.GSkyDir()`
            Sky direction of model
        name : str
            Model name
        """
        # Set point source model
        model = self._set_ptsrc(pos)

        # Set model name
        model.name(name)

        # Append model to container
        self._models.append(model)

        # Return
        return

    def _set_bkg(self, modeltype):
        """
        Set background model

        Parameters
        ----------
        modeltype : str
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
        if modeltype == 'IRF':
            model = gammalib.GCTAModelIrfBackground(spectral)
        elif modeltype == 'AEFF':
            model = gammalib.GCTAModelAeffBackground(spectral)
        elif modeltype == 'CUBE':
            model = gammalib.GCTAModelCubeBackground(spectral)
        else:
            model = None

        # Set background name
        if model is not None:
            model.name('Background')

        # Return model
        return model

    def _set_ptsrc(self, pos):
        """
        Set point source model

        Parameters
        ----------
        pos : `~gammalib.GSkyDir()`
            Sky direction of model

        Returns
        -------
        model : `~gammalib.GModelSky()`
            Point source model
        """
        # Set spatial component
        spatial = gammalib.GModelSpatialPointSource(pos)

        # Get fit parameters
        fit_pos   = self['fit_pos'].boolean()
        #fit_shape = self['fit_shape'].boolean()

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


    def _smooth_skymap(self, skymap):
        """
        Smooth the input sky map if a valid kernel was supplied

        Parameters
        ----------
        skymap : `~gammalib.GSkyMap()`
            Sky map
        """
        # Log header
        self._log_header3(gammalib.NORMAL, 'Smoothing Skymap')

        # Make sure the smoothing kernel is not 'NONE'
        if self['smoothkernel'].string().upper() != 'NONE':
            self._log_value(gammalib.NORMAL, 'Kernel',
                            self['smoothkernel'].string())
            self._log_value(gammalib.NORMAL, 'Parameter',
                            self['smoothparam'].real())
            skymap.smooth(self['smoothkernel'].string(),
                          self['smoothparam'].real())
        else:
            self._log_string(gammalib.NORMAL, 
                        'Smoothing kernel is NONE, smoothing will be ignored.')

        return


    def _load_skymap(self):
        """
        Load sky map

        Returns
        -------
        skymap : `~gammalib.GSkyMap()`
            Sky map
        """
        # Get skymap filename
        inmap = self['inmap'].filename()

        # Open sky map file
        fits = gammalib.GFits(inmap)

        # Extract primary extension as sky map
        skymap = gammalib.GSkyMap(fits.image(0))

        # Close sky map file
        fits.close()

        # Return
        return skymap.extract(0)


    def _cache_map_dirs(self, skymap):
        """
        Cache map pixel positions to save some computation time
        """
        # Setup the skymap directions
        if not skymap.is_empty():
            self._map_dirs = [skymap.inx2dir(i) for i in range(skymap.npix())]
        else:
            self._map_dirs = []


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

        # If sky map is empty then load it from input parameter
        if self._map.is_empty():
            self._map = self._load_skymap()

        self._cache_map_dirs(self._map)

        # Smooth the map
        self._smooth_skymap(self._map)

        # Write header
        self._log_header1(gammalib.NORMAL, 'Source detection')

        # Detect sources
        self._detect_sources(self._map)

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

    def models(self):
        """
        Return model container
        """
        return self._models


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cssrcdetect(sys.argv)

    # Execute application
    app.execute()
