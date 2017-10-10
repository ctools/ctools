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


def atan2(y, x):
    """
    helper function to solve triangulation problem
    """
    if x > 0:
        val = math.atan(y / x)
    elif x < 0 and y >= 0:
        val = math.atan(y / x) + math.pi
    elif x < 0 and y < 0:
        val = math.atan(y / x) - math.pi
    elif x == 0 and y > 0:
        val = math.pi / 2
    elif x == 0 and y < 0:
        val = math.pi / 2
    return val


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
        # If there are no observations in container then query the inobs
        # parameter
        if self._obs.size() == 0:
            self['inobs'].filename()

        # Query energy binning parameters
        self['ebinalg'].string()
        if self['ebinalg'].string() == 'FILE':
            self['ebinfile'].filename()
        else:
            self['emin'].real()
            self['emax'].real()
            self['enumbins'].integer()

        # Initialise background estimation method
        self._bkgmethod = self["bkgmethod"].string()

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

        # Stacking
        self._stack = self['stack'].boolean()

        # Query ahead output parameters
        if (self._read_ahead()):
            self['outroot'].string()

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
        Calculate list of reflected regions for a single observation
        :param obs: observation
        :return: list of reflected regions
        """
        pnt_dir = obs.pointing().dir()
        offset = pnt_dir.dist_deg(self._src_dir)
        print("source", self._src_dir)
        print("pointing", pnt_dir)
        print("offset", offset)
        outregions = []

        if self._srcshape == "CIRCLE":
            # angular separation of reflected regions wrt camera center
            # and number
            alpha = 2 * math.asin(self._rad / offset)
            N = int(2 * math.pi / alpha)
            alpha = 2 * math.pi / N
            # rotation angle of reference frame to generate regions
            print("x,y", self._src_dir.ra_deg() - pnt_dir.ra_deg(),
                  self._src_dir.dec_deg() - pnt_dir.dec_deg())
            theta = atan2(self._src_dir.dec_deg() - pnt_dir.dec_deg(),
                          self._src_dir.ra_deg() - pnt_dir.ra_deg())
            # theta = 2 * math.asin(math.sqrt(
            #     math.sin(self._src_dir.dec() - pnt_dir.dec()) ** 2 + math.cos(
            #         self._src_dir.dec()) * math.cos(pnt_dir.dec()) * math.sin(
            #         self._src_dir.ra() - pnt_dir.ra()) ** 2))
            print("theta", theta)
            # loop to create reflected regions
            for s in range(0, N - 1):
                print("=====")
                phi = s * alpha
                # print("phi",phi)
                xii = offset * math.cos(phi)
                yii = offset * math.sin(phi)
                print(s, xii, yii)
                xi = xii * math.cos(theta) - yii * math.sin(theta)
                yi = xii * math.sin(theta) + yii * math.cos(theta)
                print("rot", s, xi, yi)
                x = xi + pnt_dir.ra_deg()
                y = yi + pnt_dir.dec_deg()
                ctr_dir = gammalib.GSkyDir()
                ctr_dir.radec_deg(x, y)
                print(s, x, y)
                print(pnt_dir.dist_deg(ctr_dir))
                # region creation
                region = gammalib.GSkyRegionCircle(ctr_dir, self._rad)
                outregions.append(region)

        return outregions

    # Public methods
    def run(self):
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self._obs, 'Observation')

        # Set energy binning
        if self['ebinalg'] == "FILE":
            pass
            ## need to load binning from file
        else:
            etrue = gammalib.GEbounds(self["enumbins"],
                                      gammalib.GEnergy(self["emin"], 'TeV'),
                                      gammalib.GEnergy(self["emax"], 'TeV'))
            ereco = gammalib.GEbounds(self["enumbins"],
                                      gammalib.GEnergy(self["emin"], 'TeV'),
                                      gammalib.GEnergy(self["emax"], 'TeV'))

        # Loop through observations and generate pha, arf, rmf files
        outobs = gammalib.GObservations()
        for obs in self._obs:
            bkg_reg = gammalib.GSkyRegions()
            if self._bkgmethod == "REFLECTED":
                regions = self._reflected_regions(obs)
            for region in regions:
                bkg_reg.append(region)
            onoff = gammalib.GCTAOnOffObservation(obs, etrue, ereco,
                                                  self._src_region, bkg_reg)
            outobs.append(onoff)
