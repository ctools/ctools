#! /usr/bin/env python
# ==========================================================================
# Detect source in TS map
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
import sys
import glob
import os
import gammalib
import ctools
from cscripts import modutils


# ================== #
# comsrcdetect class #
# ================== #
class comsrcdetect(ctools.cscript):
    """
    Perform SRCLIX model fitting of COMPTEL observations
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._models = gammalib.GModels()
        self._fits   = gammalib.GFits()

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query parameters
        self['inmap'].filename()
        self['threshold'].real()
        self['maxsrcs'].integer()
        self['prefix'].string()

        # Query ahead output model filename
        if self._read_ahead():
            self['outmodel'].filename()
            self['outmap'].filename()
            self['outds9file'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _add_source(self, dir):
        """
        Add point source to model container

        Parameters
        ----------
        dir : `~gammalib.GSkyDir`
            Source position
        """
        # Set spatial component
        spatial = gammalib.GModelSpatialPointSource(dir)
        spatial['RA'].free()
        spatial['DEC'].free()

        # Set spectral component
        spectral = gammalib.GModelSpectralPlaw(2.0e-3, -2.0, gammalib.GEnergy(1.0,'MeV'))
        spectral['Prefactor'].min(1.0e-25)

        # Set source name
        name = '%s%3.3d' % (self['prefix'].string(), self._models.size()+1)

        # Set source model
        source = gammalib.GModelSky(spatial, spectral)
        source.name(name)
        source.tscalc(True)

        # Append source to model
        self._models.append(source)

        # Return
        return

    def _remove_pixel(self, tsmap, premap, index, threshold):
        """
        Remove pixel and all neighboring pixels above threshold from TS map

        Parameters
        ----------
        tsmap : `~gammalib.GSkyMap`
            TS map
        premap : `~gammalib.GSkyMap`
            Prefactor map
        index : int
            Pixel to be removed
        threshold : float
            TS threshold above which pixels are to be removed
        """
        # Set size of TS map
        nx = tsmap.nx()
        ny = tsmap.ny()

        # Set pixel of interest to zero
        tsmap[index] = 0.0

        # Get X and Y index of pixel
        ix0 = int(index % nx)
        iy0 = int(index / nx)

        # Remove all neighbouring pixels
        for iy in range(iy0,-1,-1):
            yoffset = iy * nx
            n       = 0
            for ix in range(ix0,-1,-1):
                if ix == ix0 and iy == iy0:
                    continue
                i = ix + yoffset
                if tsmap[i] > threshold and premap[i] > 0.0:
                    tsmap[i] = 0.0
                    n       += 1
                else:
                    break
            for ix in range(ix0+1,nx):
                if ix == ix0 and iy == iy0:
                    continue
                i = ix + yoffset
                if tsmap[i] > threshold and premap[i] > 0.0:
                    tsmap[i] = 0.0
                    n       += 1
                else:
                    break
            if n == 0:
                break
        for iy in range(iy0+1,ny):
            yoffset = iy * nx
            n       = 0
            for ix in range(ix0,-1,-1):
                if ix == ix0 and iy == iy0:
                    continue
                i = ix + yoffset
                if tsmap[i] > threshold and premap[i] > 0.0:
                    tsmap[i] = 0.0
                    n       += 1
                else:
                    break
            for ix in range(ix0+1,nx):
                if ix == ix0 and iy == iy0:
                    continue
                i = ix + yoffset
                if tsmap[i] > threshold and premap[i] > 0.0:
                    tsmap[i] = 0.0
                    n       += 1
                else:
                    break
            if n == 0:
                break

        # Return
        return

    def _search_sources(self, tsmap, premap):
        """
        Search for sources

        Parameters
        ----------
        tsmap : `~gammalib.GSkyMap`
            TS map
        premap : `~gammalib.GSkyMap`
            Prefactor map
        """
        # Write header
        self._log_header1(gammalib.NORMAL, 'Search for sources')

        # Initialise model container and diagnostics FITS file
        self._models = gammalib.GModels()
        self._fits   = gammalib.GFits()

        # Get threshold and maximum number of sources
        threshold = self['threshold'].real()
        maxsrcs   = self['maxsrcs'].integer()

        # Iterate over source finding
        while True:

            # Find pixel with largest TS value and positive Prefactor
            tsmax = 0.0
            imax  = -1
            for i in range(tsmap.npix()):
                if tsmap[i] > tsmax and premap[i] > 0.0:
                    tsmax = tsmap[i]
                    imax  = i

            # If pixel is above threshold then add point source at pixel
            if tsmax > threshold:

                # Add point source to model
                dir = tsmap.inx2dir(imax)
                self._add_source(dir)

                # Log point source
                source = 'Source %d' % (self._models.size())
                value  = '%s (TS=%.3f)' % (str(dir), tsmax)
                self._log_value(gammalib.NORMAL, source, value)

                # If maximum number of source is reached then break
                if self._models.size() >= maxsrcs:
                    break

                # Reset map pixel and all its neighbours
                self._remove_pixel(tsmap, premap, imax, threshold)

                # Save TS map after removing pixels for diagnostics
                extname = 'TS map %3.3d' % (self._models.size())
                tsmap.write(self._fits, extname)

            # ... otherwise break iterative loop
            else:
                break

        # Return
        return


    # Public methods
    def process(self):
        """
        Process the script
        """
        # Get parameters
        self._get_parameters()

        # Log header
        self._log_header1(gammalib.NORMAL, 'Input TS map')

        # Load TS and Prefactor maps
        fits   = gammalib.GFits(self['inmap'].filename())
        tsmap  = gammalib.GSkyMap(fits[0])
        premap = gammalib.GSkyMap(fits[1])

        # Search for sources
        self._search_sources(tsmap, premap)

        # Return
        return

    def save(self):
        """ 
        Save observation definition file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save source model')

        # Get output filenames
        outmodel   = self['outmodel'].filename()
        outmap     = self['outmap'].filename()
        outds9file = self['outds9file'].filename()

        # If file exists and clobber flag is false then raise an exception
        if outmodel.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outmodel.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)
        elif outmap.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outmap.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)
        elif outds9file.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outds9file.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:

            # Log model definition filename
            self._log_value(gammalib.NORMAL, 'Model definition XML file', outmodel.url())

            # Save models
            self._models.save(outmodel)

            # Log model definition filename
            self._log_value(gammalib.NORMAL, 'TS map FITS file', outmap.url())

            # Save FITS file
            self._fits.saveto(outmap, self['clobber'].boolean())

            # Log DS9 region filename
            self._log_value(gammalib.NORMAL, 'DS9 region file', outds9file.url())

            # Save DS9 region file
            modutils.models2ds9file(self._models, outds9file.url())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comsrcdetect(sys.argv)

    # Execute application
    app.execute()
