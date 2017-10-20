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
class csphagen(ctools.cscript):
    """
    Generate PHA, ARF and RMF files for classical IACT spectral analysis
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name, version
        self._name    = 'csphagen'
        self._version = ctools.__version__

        # Initialise observation container from constructor arguments
        self._obs, argv = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Initialise other variables
        self._outobs        = gammalib.GObservations()
        self._ebounds       = gammalib.GEbounds()
        self._src_dir       = gammalib.GSkyDir()
        self._src_reg       = gammalib.GSkyRegions()
        self._bkg_regs      = []
        self._excl_reg      = None
        self._has_exclusion = False
        self._srcshape      = ''
        self._rad           = 0.0

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup observations
        """
        # Setup observations (require response and allow event list, don't
        # allow counts cube)
        self._setup_observations(self._obs, True, True, False)

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
        self._srcshape = self['srcshape'].string()
        if self._srcshape == 'CIRCLE':
            self._rad = self['rad'].real()
            self._src_reg.append(gammalib.GSkyRegionCircle(self._src_dir, self._rad))

        # Query background estimation method and parameters
        bkgmethod = self['bkgmethod'].string()
        if bkgmethod == 'REFLECTED':
            self['bkgregmin'].integer()
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
        if self._obs.size() == 0:
            self._obs = self._get_observations(False)

        # Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _reflected_regions(self, obs):
        """
        Calculate list of reflected regions for a single observation (pointing)

        Parameters
        ----------
        obs : `~gammalib.GObservation()`
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
                alpha = 1.05 * 2 * self._rad / offset
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
        self._log_observations(gammalib.NORMAL, self._obs, 'Observation')

        # Write header
        self._log_header1(gammalib.NORMAL,
                          'Generation of source and background spectra')

        # Initialise run variables
        self._outobs   = gammalib.GObservations()
        self._bkg_regs = []

        # Loop through observations and generate pha, arf, rmf files
        for obs in self._obs:

            # Initialise background regions for this observation
            bkg_reg = gammalib.GSkyRegions()

            # If reflected background is requested then created reflected
            # background regions
            if self['bkgmethod'].string() == 'REFLECTED':
                bkg_reg = self._reflected_regions(obs)

            # If there are reflected regions then create On/Off observation
            # and append it to the output container
            if bkg_reg.size() >= self['bkgregmin'].integer():
                onoff = gammalib.GCTAOnOffObservation(obs, self._src_dir,
                                                      self._ebounds,
                                                      self._ebounds,
                                                      self._src_reg,
                                                      bkg_reg)
                onoff.id(obs.id())
                self._outobs.append(onoff)
                self._bkg_regs.append({'regions': bkg_reg, 'id': obs.id()})
            else:
                self._log_string(gammalib.NORMAL, 'Observation %s not included '
                                 'in spectra generation' % (obs.id()))

        # Stack observations
        if self._outobs.size() > 1 and self['stack'].boolean() == True:

            # Write header
            self._log_header1(gammalib.NORMAL, 'Stacking %d observations' %
                              (self._outobs.size()))

            # Stack observations
            stacked_obs = gammalib.GCTAOnOffObservation(self._outobs)

            # Put stacked observations in output container
            self._outobs = gammalib.GObservations()
            self._outobs.append(stacked_obs)

        # Return
        return

    def save(self):
        """ 
        Save data
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save data')

        # Get XMK output filename, prefix and clobber
        outobs  = self['outobs'].string()
        prefix  = self['prefix'].string()
        clobber = self['clobber'].boolean()

        # Loop over all observation in container
        for obs in self._outobs:

            # Set filenames
            if self['stack'].boolean():
                onname  = prefix + '_stacked_pha_on.fits'
                offname = prefix + '_stacked_pha_off.fits'
                arfname = prefix + '_stacked_arf.fits'
                rmfname = prefix + '_stacked_rmf.fits'
            elif self._outobs.size() > 1:
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
        self._outobs.save(outobs)

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

    def obs(self):
        """
        Return observation container
        """
        return self._outobs


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csphagen(sys.argv)

    # Execute application
    app.execute()
