#! /usr/bin/env python
# ==========================================================================
# Generate model for binned COMPTEL observations
#
# Copyright (C) 2019-2021 Juergen Knoedlseder
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


# ================= #
# comobsmodel class #
# ================= #
class comobsmodel(ctools.csobservation):
    """
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._models = gammalib.GModels()

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

        # Query source parameters is Right Ascension and Declination are valid
        if self['ra'].is_valid() and self['ra'].is_valid():
            self['ra'].real()
            self['dec'].real()
            self['srcname'].string()

        # Query addition components
        self['brems'].string()
        self['ic'].string()
        self['iso'].string()
        self['diffusetype'].string()
        self['bkgtype'].string()
        if self['ebinfile'].is_valid():
            self['ebinfile'].filename()
        self['bremsmap'].filename()
        self['bremscube'].filename()
        self['icmap'].filename()
        self['iccube'].filename()

        # Query ahead output model filename
        if self._read_ahead():
            self['outmodel'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _get_ebounds(self):
        """
        Get energy boundaries of observations
        """
        # Initialise energy boundaries
        ebounds = gammalib.GEbounds()

        # Loop over all observations
        for obs in self.obs():

            # Get energy range
            emin = obs.events().dre().ebounds().emin()
            emax = obs.events().dre().ebounds().emax()

            # Check if energy range is already in energy boundaries
            contained = False
            for i in range(ebounds.size()):
                if ebounds.emin(i) == emin and ebounds.emax(i) == emax:
                    contained = True
                    break

            # Append energy range if it is not yet contained in energy
            # boundaries
            if not contained:
                ebounds.append(emin, emax)

        # Return energy boundaries
        return ebounds

    def _generate_background_nodes_model(self, obs):
        """
        Generate nodes background model component for observation

        Parameters
        ----------
        obs : `~gammalib.GCOMObservation`
            COMPTEL observation

        Returns
        -------
        model : `~gammalib.GCOMModelDRBFitting`
            Background model
        """
        # Set model attributes
        name = 'Background_%s' % obs.id()
        id   = obs.id()

        # Extract attributes from observation
        drb     = obs.drb()
        npix    = drb.nchi() * drb.npsi()
        nphibar = drb.nphibar()
        dphibar = 2.0         # Fixed so far as we have no method to recover

        # Create XML model list object
        xml = gammalib.GXml()
        lib = xml.append('source_library title="source library"')
        src = lib.append('source name="%s" type="DRBFitting" instrument="COM" id="%s"' % (name, id))

        # Loop over all phibar layers
        for k in range(nphibar):

            # Compute DRB content for phibar layer
            sum_drb = 0.0
            for i in range(npix):
                index    = i + k * npix
                sum_drb += drb[index]

            # Add node if DRB content is positive
            if sum_drb > 0.0:
                phibar = (k+0.5) * dphibar
                node   = src.append('node')
                node.append('parameter name="Phibar" value="%.2f" scale="1" '
                            'min="0" max="50" free="0"' % phibar)
                node.append('parameter name="Normalization" value="1.00" '
                            'scale="1" min="0" max="1000" free="1"')

        # Create background model from XML source element
        model = gammalib.GCOMModelDRBFitting(src)

        # Return model
        return model

    def _generate_background_bins_model(self, obs):
        """
        Generate bins background model component for observation

        Parameters
        ----------
        obs : `~gammalib.GCOMObservation`
            COMPTEL observation

        Returns
        -------
        model : `~gammalib.GCOMModelDRBPhibarBins`
            Background model
        """
        # Set model attributes
        name = 'Background_%s' % obs.id()
        id   = obs.id()

        # Extract attributes from observation
        drb     = obs.drb()
        npix    = drb.nchi() * drb.npsi()
        nphibar = drb.nphibar()
        dphibar = 2.0         # Fixed so far as we have no method to recover

        # Create XML model list object
        xml = gammalib.GXml()
        lib = xml.append('source_library title="source library"')
        src = lib.append('source name="%s" type="DRBPhibarBins" instrument="COM" id="%s"' % (name, id))

        # Loop over all phibar layers
        for k in range(nphibar):

            # Compute DRB content for phibar layer
            sum_drb = 0.0
            for i in range(npix):
                index    = i + k * npix
                sum_drb += drb[index]

            # Free or fix parameter
            if sum_drb > 0.0:
                free = '1'
            else:
                free = '0'

            # Add bin
            src.append('parameter name="Normalization" value="1.00" '
                       'scale="1" min="0" max="1000" free="%s"' % free)

        # Create background model from XML source element
        model = gammalib.GCOMModelDRBPhibarBins(src)

        # Return model
        return model

    def _generate_source_model(self):
        """
        Generate source model
        """
        # Set source position
        ra        = self['ra'].real()
        dec       = self['dec'].real()
        spat      = 'PTSRC'
        spec      = 'PLAW'
        prefactor = 2.0e-3
        index     = -2.0
        name      = self['srcname'].string()

        # Set spectral component
        if spec == 'PLAW':
            spectral = gammalib.GModelSpectralPlaw(prefactor, index,
                                                   gammalib.GEnergy(1.0,'MeV'))
            spectral['Prefactor'].min(1.0e-25)

        # Set spatial component
        dir = gammalib.GSkyDir()
        dir.radec_deg(ra, dec)
        if spat == 'PTSRC':
            spatial = gammalib.GModelSpatialPointSource(dir)
        spatial['RA'].free()
        spatial['DEC'].free()

        # Set source model
        source = gammalib.GModelSky(spatial, spectral)
        source.name(name)
        source.tscalc(True)

        # Return source model
        return source

    def _recast_model(self, model, ebounds, index, min, max, fix=False):
        """
        Re-cast model

        Parameters
        ----------
        model : `~gammalib.GModelSpectral`
            Input spectral model
        ebounds : `~gammalib.GEbounds`
            Energy boundaries of observations (used if ebinfile='NONE')
        index : float
            Spectral index for bin function
        min : float
            Minimum intensities
        max : float
            Maximum intensities
        fix : boolean
            Fix intensity values

        Returns
        -------
        spectral : `~gammalib.GModelSpectralNodes` or `~gammalib.GModelSpectralBins`
            Output spectral model
        """
        # Get energy bins for the model
        if self['ebinfile'].is_valid():
            ebinfile = self['ebinfile'].filename()
            ebds     = gammalib.GEbounds(ebinfile)
        else:
            ebds = ebounds

        # Re-sample model for 'NODES' model type
        if self['diffusetype'].string() == 'NODES':

            # Setup energies vector for bin mean energies
            energies = gammalib.GEnergies()
            for i in range(ebds.size()):
                energies.append(ebds.elogmean(i))

            # Re-cast spectral model for requested energies
            spectral = gammalib.GModelSpectralNodes(model, energies)

        # Re-sample model for 'BINS' model type
        else:
            spectral = gammalib.GModelSpectralBins(model, ebds, index)

        # Set minimum and maximum values
        for i in range(ebds.size()):
            parname = 'Intensity%d' % i
            spectral[parname].min(min)
            spectral[parname].max(max)
            if fix:
                spectral[parname].fix()

        # Return spectral model
        return spectral

    def _generate_brems_model(self, ebounds):
        """
        Generate Bremsstrahlung model

        Parameters
        ----------
        ebounds : `~gammalib.GEbounds`
            Energy boundaries of observations (used if ebinfile='NONE')

        Returns
        -------
        source : `~gammalib.GModelSky`
            Bremsstrahlung model
        """
        # Append COMPTEL map if MAP is selected ...
        if self['brems'].string() == 'MAP':

            # Set spatial model (do not normalise the map)
            spatial  = gammalib.GModelSpatialDiffuseMap(self['bremsmap'].filename(), 1.0, False)

            # Set spectral model using initial scaling factors that correspond
            # to the results of Bloemen et al. (1999)
            spectral = gammalib.GModelSpectralNodes()
            spectral.append(gammalib.GEnergy(0.87, 'MeV'),  0.93)
            spectral.append(gammalib.GEnergy(1.73, 'MeV'),  0.20)
            spectral.append(gammalib.GEnergy(5.48, 'MeV'),  0.018)
            spectral.append(gammalib.GEnergy(17.32, 'MeV'), 0.0022)

            # Recast model
            spectral = self._recast_model(spectral, ebounds, -2.0, 1.0e-10, 1.0)

        # ... otherwise append Fermi-LAT GALPROP cube
        else:

            # Set spatial model
            spatial  = gammalib.GModelSpatialDiffuseCube(self['bremscube'].filename())

            # Set spectral model using initial scaling factors that correspond
            # to the results of Bloemen et al. (1999)
            spectral = gammalib.GModelSpectralNodes()
            spectral.append(gammalib.GEnergy(0.87, 'MeV'), 32.655)
            spectral.append(gammalib.GEnergy(1.73, 'MeV'), 14.232)
            spectral.append(gammalib.GEnergy(5.48, 'MeV'),  4.235)
            spectral.append(gammalib.GEnergy(17.32, 'MeV'), 1.735)

            # Recast model
            spectral = self._recast_model(spectral, ebounds, 0.0, 1.0e-5, 1.0e+5)

        # Autoscale spectral coefficients
        spectral.autoscale()

        # Set source model
        source = gammalib.GModelSky(spatial, spectral)
        source.name('Bremsstrahlung')
        source.tscalc(True)

        # Return source model
        return source

    def _generate_ic_model(self, ebounds):
        """
        Generate Inverse Compton model

        Parameters
        ----------
        ebounds : `~gammalib.GEbounds`
            Energy boundaries of observations (used if ebinfile='NONE')

        Returns
        -------
        source : `~gammalib.GModelSky`
            Inverse Compton model
        """
        # Append COMPTEL map if MAP is selected ...
        if self['ic'].string() == 'MAP':

            # Set spatial model (do not normalise the map)
            spatial  = gammalib.GModelSpatialDiffuseMap(self['icmap'].filename(), 1.0, False)

            # Set spectral model using initial scaling factors that correspond
            # to the results of Bloemen et al. (1999)
            spectral = gammalib.GModelSpectralNodes()
            spectral.append(gammalib.GEnergy(0.87, 'MeV'),  2.23)
            spectral.append(gammalib.GEnergy(1.73, 'MeV'),  0.64)
            spectral.append(gammalib.GEnergy(5.48, 'MeV'),  0.079)
            spectral.append(gammalib.GEnergy(17.32, 'MeV'), 0.011)

            # Recast model
            spectral = self._recast_model(spectral, ebounds, -2.0, 1.0e-10, 1.0)

        # ... otherwise append Fermi-LAT GALPROP cube
        else:

            # Set spatial model
            spatial  = gammalib.GModelSpatialDiffuseCube(self['iccube'].filename())

            # Set spectral model using initial scaling factors that correspond
            # to the results of Bloemen et al. (1999)
            spectral = gammalib.GModelSpectralNodes()
            spectral.append(gammalib.GEnergy(0.87, 'MeV'),  1.042)
            spectral.append(gammalib.GEnergy(1.73, 'MeV'),  0.990)
            spectral.append(gammalib.GEnergy(5.48, 'MeV'),  0.998)
            spectral.append(gammalib.GEnergy(17.32, 'MeV'), 1.239)

            # Recast model
            spectral = self._recast_model(spectral, ebounds, 0.0, 1.0e-5, 1.0e+5)

        # Autoscale spectral coefficients
        spectral.autoscale()

        # Set source model
        source = gammalib.GModelSky(spatial, spectral)
        source.name('Inverse Compton')
        source.tscalc(True)

        # Return source model
        return source

    def _generate_iso_model(self, ebounds):
        """
        Generate Isotropic model

        Parameters
        ----------
        ebounds : `~gammalib.GEbounds`
            Energy boundaries of observations (used if ebinfile='NONE')

        Returns
        -------
        source : `~gammalib.GModelSky`
            Isotropic model
        """
        # Append model if CONST is selected ...
        if 'CONST' in self['iso'].string():

            # Set spatial model, scaled so that the spectral part gives flux per
            # steradian
            spatial  = gammalib.GModelSpatialDiffuseConst(gammalib.fourpi)

            # Set spectral model using initial scaling factors that correspond
            # to Georg Weidenspointers extragalactic spectrum. Values integrated
            # over the standard energy bands communicated by Werner Collmar and
            # divided by width of standard energy bands.
            # See his e-mail from 14 July 2021.
            spectral = gammalib.GModelSpectralNodes()
            spectral.append(gammalib.GEnergy(0.87, 'MeV'),  1.327e-3/0.25)
            spectral.append(gammalib.GEnergy(1.73, 'MeV'),  2.358e-3/2.0)
            spectral.append(gammalib.GEnergy(5.48, 'MeV'),  0.658e-3/7.0)
            spectral.append(gammalib.GEnergy(17.32, 'MeV'), 0.148e-3/20.0)

            # Set fix
            fix = (self['iso'].string() == 'CONSTFIX')

            # Recast model
            spectral = self._recast_model(spectral, ebounds, -2.0, 1.0e-20, 1.0, fix=fix)

        # Autoscale spectral coefficients
        spectral.autoscale()

        # Set source model
        source = gammalib.GModelSky(spatial, spectral)
        source.name('Isotropic')
        source.tscalc(True)

        # Return source model
        return source


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
        self._log_header1(gammalib.NORMAL, 'Input observations')

        # Log input observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

        # Write header
        self._log_header1(gammalib.NORMAL, 'Energy bins')

        # Get energy boundaries
        ebounds = self._get_ebounds()

        # Log energy boundaries
        self._log_string(gammalib.NORMAL, str(ebounds))

        # Write header
        self._log_header1(gammalib.NORMAL, 'Create model components')

        # Clear output model container
        self._models.clear()

        # Append source model
        if self['ra'].is_valid() and self['ra'].is_valid():
            model = self._generate_source_model()
            self._models.append(model)

        # Optionally append Bremsstrahlung component
        if self['brems'].string() != 'NONE':
            model = self._generate_brems_model(ebounds)
            self._models.append(model)

        # Optionally append Inverse Compton component
        if self['ic'].string() != 'NONE':
            model = self._generate_ic_model(ebounds)
            self._models.append(model)

        # Optionally append Isotropic component
        if self['iso'].string() != 'NONE':
            model = self._generate_iso_model(ebounds)
            self._models.append(model)

        # Loop over all observations
        for obs in self.obs():

            # Write header
            self._log_header2(gammalib.NORMAL, self._get_obs_header(obs))

            # Log observation
            self._log_string(gammalib.NORMAL, str(obs))

            # Generate model
            if self['bkgtype'].string() == 'NODES':
                model = self._generate_background_nodes_model(obs)
            else:
                model = self._generate_background_bins_model(obs)

            # Log model
            self._log_string(gammalib.EXPLICIT, str(model))

            # Append model
            self._models.append(model)

        # Write header
        self._log_header1(gammalib.NORMAL, 'Output models')

        # Log output model
        self._log_string(gammalib.NORMAL, str(self._models))

        # Return
        return

    def save(self):
        """ 
        Save model definition XML file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save models')

        # Get output filename
        outmodel = self['outmodel'].filename()

        # If file exists and clobber flag is false then raise an exception
        if outmodel.exists() and not self['clobber'].boolean():
            msg = ('Cannot save "'+outmodel.url()+'": File already exists. '
                   'Use parameter clobber=yes to allow overwriting of files.')
            raise RuntimeError(msg)

        # Otherwise log filename and save file
        else:
            # Log filename
            self._log_value(gammalib.NORMAL, 'Model definition XML file',
                                             outmodel.url())

            # Save model definition XML file
            self._models.save(outmodel)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comobsmodel(sys.argv)

    # Execute application
    app.execute()
