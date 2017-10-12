#! /usr/bin/env python
# ==========================================================================
# Computes the PHA spectra for source/background and ARF/RMF files using the
# reflected region method
#
# Copyright (C) 2017- Luigi Tibaldo
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
    Generate pha files for classical IACT spectral analysis
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name, version
        self._name = 'csphagen'
        self._version = ctools.__version__

        # Initialise observation container from constructor arguments
        self._obs, argv = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

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

        # Initialise background estimation method
        self._bkgmethod = self["bkgmethod"].string()
        if self._bkgmethod == "REFLECTED":
            self._bkgregmin = self["bkgregmin"].integer()
        self._maxoffset = self["maxoffset"].real()

        # Initialise source position/region querying relevant parameters
        self._src_dir = gammalib.GSkyDir()
        self._src_reg = gammalib.GSkyRegions()
        coordsys = self['coordsys'].string()
        if coordsys == "CEL":
            ra = self['ra'].real()
            dec = self['dec'].real()
            self._src_dir.radec_deg(ra, dec)
        elif coordsys == "GAL":
            glon = self['glon'].real()
            glat = self['glat'].real()
            self._src_dir.lb_deg(glon, glat)
        self._srcshape = self['srcshape'].string()
        if self._srcshape == "CIRCLE":
            self._rad = self['rad'].real()
            self._src_reg.append(
                gammalib.GSkyRegionCircle(self._src_dir, self._rad))

        # exclusion map
        if self["exclusion"].filename().is_fits():
            self._excl_reg = gammalib.GSkyRegionMap(
                self["exclusion"].filename())
            self._has_exclusion = True

        # Stacking
        self._stack = self['stack'].boolean()

        # Query ahead output parameters
        if (self._read_ahead()):
            self._outroot = self['outroot'].string()

        # Set some fixed parameters
        self._chatter = self["chatter"].integer()
        self._clobber = self["clobber"].boolean()
        self._debug = self["debug"].boolean()

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
        :param obs: observation
        :return: list of reflected regions
        """
        outregions = []
        pnt_dir = obs.pointing().dir()
        offset = pnt_dir.dist_deg(self._src_dir)
        if offset <= self._rad or offset >= self["maxoffset"].real():
            if self._logExplicit():
                msg = 'Observation {} pointed at {} deg from source\n'.format(
                    obs.id(), offset)
                self._log(msg)
        else:
            posang = pnt_dir.posang_deg(self._src_dir)
            if self._srcshape == "CIRCLE":
                # angular separation of reflected regions wrt camera center
                # and number
                alpha = 1.05 * 2 * self._rad / offset
                # 1.05 ensures background regions do not overlap due to numerical precision issues
                N = int(2 * math.pi / alpha)
                if N < self._bkgregmin + 3:
                    if self._logExplicit():
                        msg = 'Observation {}: insufficient regions for background estimation\n'.format(
                            obs.id(), offset)
                        self._log = (msg)
                else:
                    alpha = 360. / N
                    # loop to create reflected regions
                    for s in range(2, N - 1):
                        dphi = s * alpha
                        ctr_dir = pnt_dir.clone()
                        ctr_dir.rotate_deg(posang + dphi, offset)
                        region = gammalib.GSkyRegionCircle(ctr_dir, self._rad)
                        if self._has_exclusion:
                            if self._excl_reg.overlaps(region):
                                if self._logVerbose():
                                    msg = 'Observation {}: reflected region overlaps with exclusion region\n'.format(
                                        obs.id(), offset)
                                    self._log(msg)
                            else:
                                outregions.append(region)
                        else:
                            outregions.append(region)

        return outregions

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

        self._log_header1(gammalib.NORMAL,
                          'Generation of source and background spectra')
        # Loop through observations and generate pha, arf, rmf files
        outobs = gammalib.GObservations()
        etrue = self._ebounds
        ereco = self._ebounds
        for obs in self._obs:
            bkg_reg = gammalib.GSkyRegions()
            if self._bkgmethod == "REFLECTED":
                regions = self._reflected_regions(obs)
            for region in regions:
                bkg_reg.append(region)
            if bkg_reg.size() >= 1 or (
                            self._bkgmethod == "REFLECTED" and bkg_reg.size() >= self._bkgregmin):
                onoff = gammalib.GCTAOnOffObservation(obs, etrue, ereco,
                                                      self._src_reg, bkg_reg)
                onoff.id(obs.id())
                outobs.append(onoff)
            else:
                if self._logNormal():
                    msg = 'Observation {} not included in spectra generation\n'.format(
                        obs.id())
                    self._log(msg)

        # Save PHA, ARF and RMFs
        for obs in outobs:
            obs.on_spec().save(
                self._outroot + '_{}_pha_on.fits'.format(obs.id()),
                True)
            obs.off_spec().save(
                self._outroot + '_{}_pha_off.fits'.format(obs.id()),
                True)
            obs.arf().save(self._outroot + '_{}_arf.fits'.format(obs.id()),
                           True)
            obs.rmf().save(self._outroot + '_{}_rmf.fits'.format(obs.id()),
                           True)
            obs.on_regions().save(self._outroot + '_{}_on.reg'.format(obs.id()))
            obs.off_regions().save(
                self._outroot + '_{}_off.reg'.format(obs.id()))

        # Save On/Off observations
        outname = self._outroot + '.xml'
        outobs.save(outname)
        # Log filename
        if self._logNormal():
            self._log(
                'Output observation definition XML file: {}\n'.format(outname))

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

        # Save residual map
        self.save()

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
