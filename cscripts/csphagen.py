#! /usr/bin/env python
# ==========================================================================
# Computes the PHA spectra for source/background and ARF/RMF files using the
# reflected region method
#
# Copyright (C) 2017 Luigi Tibaldo
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
import gammalib
import ctools
import math
import sys


# =============== #
# csfindobs class #
# =============== #
class csphagen(ctools.csobservation):
    """
    Generate PHA, ARF and RMF files for classical IACT spectral analysis
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise other variables
        self._ebounds       = gammalib.GEbounds()
        self._src_dir       = gammalib.GSkyDir()
        self._src_reg       = gammalib.GSkyRegions()
        self._bkg_regs      = []
        self._excl_reg      = None
        self._has_exclusion = False
        self._srcshape      = ''
        self._rad           = 0.0
        self._bkg_reg_files = {}     # dictionary: bkg region file per observation {obsid:filepath, ..}

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup observations

        Raises
        ------
        RuntimeError
            Multiple On-Regions in on region file or missing Off-Regions files.
        """
        # Setup observations (require response and allow event list, don't
        # allow counts cube)
        self._setup_observations(self.obs(), True, True, False)

        # Set energy bounds
        self._ebounds = self._create_ebounds()

        # Initialise source position/region querying relevant parameters
        self._src_dir = gammalib.GSkyDir()
        self._src_reg = gammalib.GSkyRegions()
        coordsys = self['coordsys'].string()
        if coordsys == 'CEL':
            ra  = self['ra'].real()
            dec = self['dec'].real()
            self._src_dir.radec_deg(ra, dec)
        elif coordsys == 'GAL':
            glon = self['glon'].real()
            glat = self['glat'].real()
            self._src_dir.lb_deg(glon, glat)

        # Query background estimation method and parameters
        bkgmethod = self['bkgmethod'].string()
        if bkgmethod == 'REFLECTED':

            self._srcshape = self['srcshape'].string()
            if self._srcshape == 'CIRCLE':
                self._rad = self['rad'].real()
                self._src_reg.append(gammalib.GSkyRegionCircle(self._src_dir, self._rad))

            self['bkgregmin'].integer()

        elif bkgmethod == 'CUSTOM':

            obs_file     = self['inobs'].filename()
            src_reg_file = gammalib.GFilename()

            # obs def XML: look for region files
            if obs_file.file()[-3:] == 'xml':
                # open and load XML obs definition file
                xml = gammalib.GXml()
                xml.load(obs_file)

                # Go through elements in <observationlist>
                for elem in xml[0]:

                    # Get observation id from <observation> header
                    obs_id = elem.attribute('id')

                    # Get region file path from <observation> body
                    off_regions_file = ''
                    on_region_file   = ''
                    for entry in elem:

                        # <parameter name="OffRegions" file=".."/> contains region file
                        if entry.attribute('name') == 'OffRegions':
                            off_regions_file = entry.attribute('file')

                        # <parameter name="OnRegion" file=".."/> contains region file.
                        if entry.attribute('name') == 'OnRegion':
                            on_region_file = entry.attribute('file')

                    # Store file paths
                    self._bkg_reg_files[obs_id] = gammalib.GFilename(off_regions_file)

                    if on_region_file != '':
                        src_reg_file = gammalib.GFilename(on_region_file)

                    # Not existing bkg region files are set to None (for every obs)
                    if not self._bkg_reg_files[obs_id].exists():
                        self._bkg_reg_files[obs_id] = None

                # Not existing src region files are set to None (one time)
                if not src_reg_file.exists():
                    src_reg_file = None

            else:
                # event FITS files: query later 1 bkg & 1 src region file
                self._bkg_reg_files = {0:None}
                src_reg_file  = None

            # Query missing bkg region file
            if len([True for x in self._bkg_reg_files if self._bkg_reg_files[x] is None]) > 0:
                if len(self._bkg_reg_files) > 1:
                    # Querying of multiple regions not supported
                    raise RuntimeError('There are missing bkg region files in the obs definition'
                                       +' XML. Please specify all bkg region files!')
                else:
                    # Query bkg region file for single observation
                    self._bkg_reg_files[self._bkg_reg_files.keys()[0]] = self['bkgreg'].filename()

            # Query missing src region file
            if src_reg_file is None:
                src_reg_file  = self['srcreg'].filename()

            # Set up src (=On) region. Only 1 region is allowed
            self._src_reg = gammalib.GSkyRegions(src_reg_file)

            if len(self._src_reg) > 1:
                raise RuntimeError('Only 1 ON region is allowed.')

        self['maxoffset'].real()

        # Exclusion map
        if self['inexclusion'].filename().is_fits():
            self._excl_reg = gammalib.GSkyRegionMap(self['inexclusion'].filename())
            self._has_exclusion = True
        else:
            self._has_exclusion = False

        # Query remaining parametersStacking
        self['stack'].boolean()

        # Query ahead output parameters
        if (self._read_ahead()):
            self['outobs'].string()
            self['prefix'].string()

        # If there are no observations in container then get them from the
        # parameter file
        if self.obs().is_empty():
            self.obs(self._get_observations(False))

        # Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _reflected_regions(self, obs):
        """
        Calculate list of reflected regions for a single observation (pointing)

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation()`
            CTA observation

        Returns
        -------
        regions : `~gammalib.GSkyRegions`
            List of reflected regions
        """
        # Initialise list of reflected regions
        regions = gammalib.GSkyRegions()

        # Get offset angle of source
        pnt_dir = obs.pointing().dir()
        offset  = pnt_dir.dist_deg(self._src_dir)

        # If ...
        if offset <= self._rad or offset >= self['maxoffset'].real():
            self._log_string(gammalib.EXPLICIT, 'Observation %s pointed at %.3f '
                             'deg from source' % ((obs.id(), offset)))

        # ... otherwise
        else:
            posang = pnt_dir.posang_deg(self._src_dir)
            if self._srcshape == 'CIRCLE':
                # angular separation of reflected regions wrt camera center
                # and number
                alpha = 1.05 * 2.0 * self._rad / offset
                # 1.05 ensures background regions do not overlap due to
                # numerical precision issues
                N = int(2.0 * math.pi / alpha)
                if N < self['bkgregmin'].integer() + 3:
                    self._log_string(gammalib.EXPLICIT, 'Observation %s: '
                                     'insufficient regions for background '
                                     'estimation' % (obs.id()))

                # Otherwise loop over position angle to create reflected
                # regions
                else:
                    alpha = 360.0 / N
                    for s in range(2, N - 1):
                        dphi    = s * alpha
                        ctr_dir = pnt_dir.clone()
                        ctr_dir.rotate_deg(posang + dphi, offset)
                        region = gammalib.GSkyRegionCircle(ctr_dir, self._rad)
                        if self._has_exclusion:
                            if self._excl_reg.overlaps(region):
                                self._log_string(gammalib.VERBOSE, 'Observation '
                                                 '%s: reflected region overlaps '
                                                 'with exclusion region' %
                                                 (obs.id()))
                            else:
                                regions.append(region)
                        else:
                            regions.append(region)

        # Return reflected regions
        return regions

    def _regions_from_file(self, obs, file_path):
        """
        Read list of reflected regions for a single observation (pointing) from file.
        Perform basic tests for circle shaped regions.

        Parameters
        ----------
        obs : `~gammalib.GCTAObservation()`
            CTA observation
        file_path: `~gammalib.GFilename()`
            File path to file containing off regions.

        Returns
        -------
        regions : `~gammalib.GSkyRegions`
            List of reflected regions
        """

        # Load regions from ds9 file of FITS WCS map
        regions  = gammalib.GSkyRegions( file_path )

        # Do region tests only for circle shaped on-region. Else simply return list.
        on_region = self._src_reg[0]
        if isinstance(on_region, gammalib.GSkyRegionCircle):

            # Prepare parameters for comparison
            pnt_dir  = obs.pointing().dir()
            offset_on = pnt_dir.dist_deg(self._src_dir)
            radius_on = on_region.radius()

            # Do some tests
            for region in regions:

                # Do region tests only for circle shaped off-regions
                if isinstance(region, gammalib.GSkyRegionCircle):

                    # Gather parameters
                    offset_off = pnt_dir.dist_deg(region.centre())
                    radius_off = region.radius()

                    # Warn: different radius
                    if radius_off != radius_on:
                        self._log_string(gammalib.EXPLICIT, 'Off region radius differs from on'
                                         'region radius: %.2f vs %.2f' % ((radius_off, radius_on)))

                    # Warn: offset is too small or too large
                    if offset_off <= radius_on or offset_off >= self['maxoffset'].real():
                        self._log_string(gammalib.EXPLICIT, 'Region offset invalid for observation '
                                         '%s (%.3f deg) ' % ((obs.id(), offset_off)))

                    # Warn: offset of off and on region differ
                    if offset_off != offset_on:
                        self._log_string(gammalib.EXPLICIT, 'Off region offset differs '
                                         '%.3f deg from on region offset' % (offset_off-offset_on))

                    # end: off-regions test for circle shaped regions
                # end: for loop iterating over off regions
            # end: do region tests only if on-region is circle shaped

        # Inform user how many regions were loaded
        self._log_string(gammalib.NORMAL, 'Read %i regions from file %s.' % (regions.size(), file_path.url()))

        # return list of regions
        return regions

    def _set_models(self, obs):
        """
        Set models in observation container

        The method replaces all "CTA" background models by "CTAOnOff"
        background models.

        Parameters
        ----------
        obs : `~gammalib.GObservations()`
            Observation container

        Returns
        -------
        obs : `~gammalib.GObservations()`
            Observation container
        """
        # Initialise model container
        models = gammalib.GModels()

        # Loop over all models and replace CTA background model by CTAOnOff
        # background model
        for model in self.obs().models():
            if 'GCTA' in model.classname() and 'CTA' in model.instruments():
                model.instruments('CTAOnOff')
            models.append(model)

        # Append model to observation container
        obs.models(models)

        # Return observation container
        return obs

    def _etrue_ebounds(self):
        """
        Set true energy boundaries

        Returns
        -------
        ebounds : `~gammalib.GEbounds()`
            True energy boundaries
        """
        # Determine minimum and maximum energies
        emin = self._ebounds.emin() * 0.5
        emax = self._ebounds.emax() * 1.2
        if emin.TeV() < self['etruemin'].real():
            emin = gammalib.GEnergy(self['etruemin'].real(), 'TeV')
        if emax.TeV() > self['etruemax'].real():
            emax = gammalib.GEnergy(self['etruemax'].real(), 'TeV')

        # Determine number of energy bins
        n_decades = (emax.log10TeV() - emin.log10TeV())
        n_bins    = int(n_decades * float(self['etruebins'].integer())) + 1

        # Set energy boundaries
        ebounds = gammalib.GEbounds(n_bins, emin, emax)

        # Write header
        self._log_header1(gammalib.TERSE, 'True energy binning')

        # Log true energy bins
        for i in range(ebounds.size()):
            value = '%s - %s' % (str(ebounds.emin(i)), str(ebounds.emax(i)))
            self._log_value(gammalib.TERSE, 'Bin %d' % (i+1), value)

        # Return energy boundaries
        return ebounds


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

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Observation')

        # Set true energy bins
        etrue_ebounds = self._etrue_ebounds()

        # Write header
        self._log_header1(gammalib.TERSE, 'Spectral binning')

        # Log reconstructed energy bins
        for i in range(self._ebounds.size()):
            value = '%s - %s' % (str(self._ebounds.emin(i)),
                                 str(self._ebounds.emax(i)))
            self._log_value(gammalib.TERSE, 'Bin %d' % (i+1), value)

        # Write header
        self._log_header1(gammalib.NORMAL,
                          'Generation of source and background spectra')

        # Initialise run variables
        outobs         = gammalib.GObservations()
        self._bkg_regs = []

        # Loop through observations and generate pha, arf, rmf files
        for obs in self.obs():

            # Log current observation
            self._log_string(gammalib.NORMAL, 'Process Observation %s' % (obs.id()))

            # Initialise background regions for this observation
            bkg_reg = gammalib.GSkyRegions()

            # If reflected background is requested then created reflected
            # background regions
            if self['bkgmethod'].string() == 'REFLECTED':
                bkg_reg = self._reflected_regions(obs)
            elif self['bkgmethod'].string() == 'CUSTOM':
                # Get file path of bkg region file from dictionary
                if len(self._bkg_reg_files) > 1:
                    bkg_reg_file = self._bkg_reg_files[obs.id()]
                else:
                    # for inobs==event list no obs id is set.
                    bkg_reg_file = self._bkg_reg_files[self._bkg_reg_files.keys()[0]]

                # Load off regions from region file for current observation
                bkg_reg = self._regions_from_file(obs, bkg_reg_file)


            # If there are reflected regions then create On/Off observation
            # and append it to the output container
            if bkg_reg.size() >= self['bkgregmin'].integer():
                onoff = gammalib.GCTAOnOffObservation(obs, self._src_dir,
                                                      etrue_ebounds,
                                                      self._ebounds,
                                                      self._src_reg,
                                                      bkg_reg)
                onoff.id(obs.id())
                outobs.append(onoff)
                self._bkg_regs.append({'regions': bkg_reg, 'id': obs.id()})
            else:
                self._log_string(gammalib.NORMAL, 'Observation %s not included '
                                 'in spectra generation' % (obs.id()))

        # Stack observations
        if outobs.size() > 1 and self['stack'].boolean() == True:

            # Write header
            self._log_header1(gammalib.NORMAL, 'Stacking %d observations' %
                              (outobs.size()))

            # Stack observations
            stacked_obs = gammalib.GCTAOnOffObservation(outobs)

            # Put stacked observations in output container
            outobs = gammalib.GObservations()
            outobs.append(stacked_obs)

        # Set models in output container
        outobs = self._set_models(outobs)

        # Set observation container
        self.obs(outobs)

        # Return
        return

    def save(self):
        """ 
        Save data
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save data')

        # Get XML output filename, prefix and clobber
        outobs  = self['outobs'].string()
        prefix  = self['prefix'].string()
        clobber = self['clobber'].boolean()

        # Loop over all observation in container
        for obs in self.obs():

            # Set filenames
            if self['stack'].boolean():
                onname  = prefix + '_stacked_pha_on.fits'
                offname = prefix + '_stacked_pha_off.fits'
                arfname = prefix + '_stacked_arf.fits'
                rmfname = prefix + '_stacked_rmf.fits'
            elif self.obs().size() > 1:
                onname  = prefix + '_%s_pha_on.fits' % (obs.id())
                offname = prefix + '_%s_pha_off.fits' % (obs.id())
                arfname = prefix + '_%s_arf.fits' % (obs.id())
                rmfname = prefix + '_%s_rmf.fits' % (obs.id())
            else:
                onname  = prefix + '_pha_on.fits'
                offname = prefix + '_pha_off.fits'
                arfname = prefix + '_arf.fits'
                rmfname = prefix + '_rmf.fits'

            # Save files
            obs.on_spec().save(onname, clobber)
            obs.off_spec().save(offname, clobber)
            obs.arf().save(arfname, clobber)
            obs.rmf().save(rmfname, clobber)

            # Log file names
            self._log_value(gammalib.NORMAL, 'PHA on file', onname)
            self._log_value(gammalib.NORMAL, 'PHA off file', offname)
            self._log_value(gammalib.NORMAL, 'ARF file', arfname)
            self._log_value(gammalib.NORMAL, 'RMF file', arfname)

        # Save observation definition XML file
        self.obs().save(outobs)

        # Log file name
        self._log_value(gammalib.NORMAL, 'Obs. definition XML file', outobs)

        # Save ds9 On region file
        regname = prefix + '_on.reg'
        self._src_reg.save(regname)
        self._log_value(gammalib.NORMAL, 'On region file', regname)

        # Save ds9 Off region files
        for bkg_reg in self._bkg_regs:

            # Set filename
            if len(self._bkg_regs) > 1:
                regname = prefix + '_%s_off.reg' % (bkg_reg['id'])
            else:
                regname = prefix + '_off.reg'

            # Save ds9 region file
            bkg_reg['regions'].save(regname)

            # Log file name
            self._log_value(gammalib.NORMAL, 'Off region file', regname)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csphagen(sys.argv)

    # Execute application
    app.execute()
