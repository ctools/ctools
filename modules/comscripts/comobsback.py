#! /usr/bin/env python
# ==========================================================================
# Generate background model for COMPTEL observations
#
# Copyright (C) 2021 Juergen Knoedlseder
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
# comobsback class #
# ================ #
class comobsback(ctools.csobservation):
    """
    Generate background model for COMPTEL observations
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
        self['outfolder'].string()
        self['bkgmethod'].string()

        # Check parameters for BGDLIXA
        if self['bkgmethod'].string() == 'BGDLIXA':

            # Get parameters
            nrunav = self['nrunav'].integer()
            navgr  = self['navgr'].integer()
            nincl  = self['nincl'].integer()
            nexcl  = self['nexcl'].integer()

            # Check for incorrect parameters
            if nexcl < 0 or nexcl >= nincl:
                msg = 'Incorrect value %d for nexcl (bins to exclude).' % nexcl
                raise RuntimeError(msg)
            if nexcl != 0 and 2*int(nexcl/2) == nexcl :
                msg = 'nexcl=%d (bins to exclude) should be zero or odd number.' % nexcl
                raise RuntimeError(msg)
            if nincl < 3 or 2*int(nincl/2) == nincl:
                msg = 'nincl=%d (bins to include) should be odd and >= 3.' % nincl
                raise RuntimeError(msg)
            if navgr < 1 or 2*int(navgr/2) == navgr :
                msg = 'navgr=%d should be odd and >= 1.' % navgr
                raise RuntimeError(msg)

        # Query ahead output model filename
        if self._read_ahead():
            self['outobs'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _generate_drb_phinor(self, drename, drgname, drm):
        """
        Generate DRB by Phibar-normalisation of DRG to DRE

        Parameters
        ----------
        drename : str
            DRE filename
        drgname : str
            DRG filename
        drm : `~gammalib.GCOMDri`
            Source DRM

        Returns
        -------
        drbname : str
            DRB filename
        """
        # Get suffix
        suffix    = self['suffix'].string()
        outfolder = self['outfolder'].string()

        # Set DRB filename
        drbname = '%s/%s' % (outfolder, os.path.basename(drename.replace('dre',
                             'drb-phinor%s' % (suffix))))
        drbfile = gammalib.GFilename(drbname)

        # Write header
        self._log_header3(gammalib.NORMAL, 'Compute DRB using Phibar normalisation')

        # If DRB file exists then do nothing
        if drbfile.exists():
             self._log_value(gammalib.NORMAL, 'DRB file exists', drbfile.url())

        # ... otherwise compute and save it
        else:

            # Set DRE and DRG filenames
            drefile = gammalib.GFilename(drename)
            drgfile = gammalib.GFilename(drgname)

            # Load DRE and DRG
            dre = gammalib.GCOMDri(drefile)
            drb = gammalib.GCOMDri(drgfile)

            # Subtract DRM
            sum_dre = 0.0
            sum_drm = 0.0
            for i in range(dre.size()):
                sum_dre += dre[i]
                sum_drm += drm[i]
                dre[i]  -= drm[i]
            self._log_value(gammalib.NORMAL, 'Events in DRE', sum_dre)
            self._log_value(gammalib.NORMAL, 'Events in DRM', sum_drm)

            # Get dataspace dimensions
            nchi    = dre.nchi()
            npsi    = dre.npsi()
            nphibar = dre.nphibar()
            npix    = nchi * npsi

            # Multiply DRG pixels by solid angle
            for i in range(npix):
                omega = drb.map().solidangle(i)
                for i2 in range(nphibar):
                    index       = i + i2 * npix
                    drb[index] *= omega

            # Phibar Normalise DRB
            for i2 in range(nphibar):
                sum_dre = 0.0
                sum_drb = 0.0
                for i in range(npix):
                    index    = i + i2*npix
                    sum_dre += dre[index]
                    sum_drb += drb[index]
                if sum_drb > 0:
                    for i in range(npix):
                        index       = i + i2*npix
                        drb[index] *= sum_dre / sum_drb

            # Save DRB
            drb.save(drbfile, self['clobber'].boolean())

            # Log creation
            self._log_value(gammalib.NORMAL, 'DRB file created', drbfile.url())

        # Return
        return drbname

    def _generate_drb_bgdlixA(self, drename, drgname, drm,
                              average_dre_by_drg=True, debug=False):
        """
        Generate DRB from DRG and DRE file with BGDLIX method A

        The method generates a background model using the BGDLIX method A,
        taking into account a source model DRM so that the method can be
        used for an iterative SRCLIX approach.

        Parameters
        ----------
        drename : str
            DRE filename
        drgname : str
            DRG filename
        drm : `~gammalib.GCOMDri`
            Source DRM
        average_dre_by_drg : bool
            Compute sum_i DRE / sum_i DRG instead of sum_i DRE/DRG (default: True)
        debug : bool
            Write some files for debugging (default: False)

        Returns
        -------
        drbname : str
            DRB filename
        """
        # Get suffix and outfolder
        suffix    = self['suffix'].string()
        outfolder = self['outfolder'].string()

        # Get task parameters
        nrunav = self['nrunav'].integer()
        navgr  = self['navgr'].integer()
        nincl  = self['nincl'].integer()
        nexcl  = self['nexcl'].integer()

        # Set DRB filename
        drbname = '%s/%s' % (outfolder, os.path.basename(drename.replace('dre',
                             'drb-bgdlixA-nr%d-na%d-ni%d-ne%d%s' %
                             (nrunav, navgr, nincl, nexcl, suffix))))
        drbfile = gammalib.GFilename(drbname)

        # Write header
        self._log_header3(gammalib.NORMAL, 'Compute DRB using BGDLIX method A')
        self._log_value(gammalib.NORMAL, 'Bins for running average', nrunav)
        self._log_value(gammalib.NORMAL, 'Bins used for averaging', navgr)
        self._log_value(gammalib.NORMAL, 'Phibar layers to include', nincl)
        self._log_value(gammalib.NORMAL, 'Phibar layers to exclude', nexcl)

        # If DRB file exists then do nothing
        if drbfile.exists():
            self._log_value(gammalib.NORMAL, 'DRB file exists', drbfile.url())

        # ... otherwise compute and save it
        else:

            # Set DRE and DRG filenames
            drefile = gammalib.GFilename(drename)
            drgfile = gammalib.GFilename(drgname)

            # Load DRE and DRG
            dre = gammalib.GCOMDri(drefile)
            drg = gammalib.GCOMDri(drgfile)

            # Get dataspace dimensions
            nchi    = dre.nchi()
            npsi    = dre.npsi()
            nphibar = dre.nphibar()
            npix    = nchi * npsi

            # Multiply DRG pixels by solid angle
            for i in range(npix):
                omega = drg.map().solidangle(i)
                for i2 in range(nphibar):
                    index       = i + i2 * npix
                    drg[index] *= omega

            # Precompute half running average lengths
            navgr2 = int(navgr/2)
            nexcl2 = int(nexcl/2)
            nincl2 = int(nincl/2)

            # Initialise SCRAT and DRB
            scrat = drg.copy()
            drb   = drg.copy()

            # Do Phibar normalization of DRG and store result in SCRAT.
            #
            # Equation (3.12) in Rob van Dijk's thesis (SCRAT = B^0_L)
            for i3 in range(nphibar):
                sum_dre = 0.0
                sum_drg = 0.0
                sum_drm = 0.0
                for i2 in range(npsi):
                    for i1 in range(nchi):
                        ipix     = i3*npix + i2*nchi + i1
                        sum_dre += dre[ipix]
                        sum_drg += drg[ipix]
                        sum_drm += drm[ipix]
                if sum_drg > 0.0:
                    scale = (sum_dre-sum_drm) / sum_drg
                    for i2 in range(npsi):
                        for i1 in range(nchi):
                            ipix        = i3*npix + i2*nchi + i1
                            scrat[ipix] = drg[ipix] * scale

            # Debug: save B^0_L
            if debug:
                srcatname = drbname.replace('drb', 'scrat1')
                srcatfile = gammalib.GFilename(srcatname)
                scrat.save(srcatfile, self['clobber'].boolean())

            # Do 3D running average of DRE and Phibar normalised DRG (SCRAT),
            # determine ratio and multiply the result by SCRAT. The result
            # is a Phibar normalised DRG that is normalised over the 3D region
            # to the DRE. The 3D region is all Phibar layers and "nrunav"
            # Chi/Psi pixels around the pixel of consideration.
            #
            # Equation (3.13) in Rob van Dijk's thesis (SCRAT = B^0C_L)
            if nrunav >= 1:
                for i2 in range(npsi):
                    for i1 in range(nchi):
                        sum_dre   = 0.0
                        sum_drm   = 0.0
                        sum_scrat = 0.0
                        for j3 in range(nphibar):
                            for j2 in range(max(0, i2 - nrunav), min(npsi, i2 + nrunav + 1)):
                                for j1 in range(max(0, i1 - nrunav), min(nchi, i1 + nrunav + 1)):
                                    ipix = j3*npix + j2*nchi + j1
                                    if drg[ipix] != 0.0:
                                        sum_dre   += dre[ipix]
                                        sum_drm   += drm[ipix]
                                        sum_scrat += scrat[ipix]
                        if sum_scrat != 0.0:
                            scale = (sum_dre-sum_drm) / sum_scrat
                            for i3 in range(nphibar):
                                ipix        = i3*npix + i2*nchi + i1
                                scrat[ipix] = scrat[ipix] * scale

            # Debug: save B^0C_L
            if debug:
                srcatname = drbname.replace('drb', 'scrat2')
                srcatfile = gammalib.GFilename(srcatname)
                scrat.save(srcatfile, self['clobber'].boolean())

            # Pre-compute the DRB by adjusting the SCRAT array over "nincl"
            # Phibar layers
            #
            # Part of Equation (3.14) in Rob van Dijk's thesis (DRB = B^0C_L / sum B^0C_L)
            for j3 in range(nphibar):

                # Determine index range for sum. There are two sums that
                # go over [isel1, iex1[ and [ixe2, isel2[
                iex1 = j3 - nexcl2
                iex2 = j3 + nexcl2 + 1
                if nexcl == 0:
                    iex2 -= 1
                isel1 = max((j3 - nincl2), 0)
                isel2 = min((j3 + nincl2 + 1), nphibar)
                if isel1 == 0:
                    isel2 = min(nincl, nphibar)
                if isel2 == nphibar:
                    isel1 = max(nphibar - nincl, 0)

                # Loop over all Chi/Psi pixels
                for i2 in range(npsi):
                    for i1 in range(nchi):

                        # Precompute DRB by dividing SCRAT by the running
                        # average over Phibar
                        ipix      = j3*npix + i2*nchi + i1
                        drb[ipix] = 0.0
                        if scrat[ipix] != 0.0:
                            sum_scrat = 0.0
                            if iex1 >= 1:
                                for i3 in range(isel1, iex1):
                                    ipixu      = i3*npix + i2*nchi + i1
                                    sum_scrat += scrat[ipixu]
                            if iex2 < nphibar:
                                for i3 in range(iex2, isel2):
                                    ipixu      = i3*npix + i2*nchi + i1
                                    sum_scrat += scrat[ipixu]
                            if sum_scrat != 0.0:
                                drb[ipix] = scrat[ipix] / sum_scrat

            # Debug: save B^0C_L / sum B^0C_L
            if debug:
                srcatname = drbname.replace('drb', 'scrat3')
                srcatfile = gammalib.GFilename(srcatname)
                drb.save(srcatfile, self['clobber'].boolean())

            # Make background model
            #
            # Part of Equation (3.14) in Rob van Dijk's thesis that corresponds
            # to G * [E / G]^s that will be stored in SCRAT.
            #
            # NOTES:
            # - we replaced DRE by DRE-DRM in the nominator of the running
            #   average since and source component should of course be
            #   subtracted
            # - we implemented an alternative running average computation where
            #   the ratio between the DRE and DRG sums is taken. This avoids
            #   the divergence of the ratio at the edge of the DRG. This
            #   algorithm is enabled by setting average_dre_by_drg=True
            scrat = dre.copy()
            if navgr != 1:
                for i3 in range(nphibar):
                    for i2 in range(npsi):
                        for i1 in range(nchi):

                            # Compute average DRE/DRG over small patch in Chi/Psi
                            if average_dre_by_drg:

                                # Compute average using sum_i DRE / sum_i DRG
                                sum_dre = 0.0
                                sum_drm = 0.0
                                sum_drg = 0.0
                                for j2 in range(max(0, i2 - navgr2), min(npsi, i2 + navgr2 + 1)):
                                    for j1 in range(max(0, i1 - navgr2), min(nchi, i1 + navgr2 + 1)):
                                        i        = i3*npix + j2*nchi + j1
                                        sum_dre += dre[i]
                                        sum_drm += drm[i]
                                        sum_drg += drg[i]

                                # Multiply average by DRG and store in SCRAT
                                ipix = i3*npix + i2*nchi + i1
                                if sum_drg > 0.0:
                                    scrat[ipix] = drg[ipix] * (sum_dre-sum_drm)/sum_drg
                                else:
                                    scrat[ipix] = 0.0

                            # ... otherwise use original algorithm
                            else:
                            
                                # Compute average sum_i DRE / DRG
                                sum1 = 0.0
                                nsum = 0
                                for j2 in range(max(0, i2 - navgr2), min(npsi, i2 + navgr2 + 1)):
                                    for j1 in range(max(0, i1 - navgr2), min(nchi, i1 + navgr2 + 1)):
                                        i = i3*npix + j2*nchi + j1
                                        if drg[i] != 0.0:
                                            sum1 += dre[i]/drg[i]
                                            nsum += 1

                                # Multiply average by DRG and store in SCRAT
                                ipix = i3*npix + i2*nchi + i1
                                if nsum != 0:
                                    scrat[ipix] = drg[ipix] * sum1/float(nsum)
                                else:
                                    scrat[ipix] = 0.0

            # Debug: save scrat
            if debug:
                srcatname = drbname.replace('drb', 'scrat4')
                srcatfile = gammalib.GFilename(srcatname)
                scrat.save(srcatfile, self['clobber'].boolean())

            # Make background model by summing SCRAT over "nincl" Phibar
            # layers and by multiplying the result with the DRB
            #
            # Last part of Equation (3.14) in Rob van Dijk's thesis that sums
            # G * [E / G]^s over a subset of Phibar layers (the first parenthesis)
            # and multplies with B^0C_L / sum B^0C_L.
            for j3 in range(nphibar):

                # Determine index range for sum. There are two sums than
                # go over [isel1, iex1[ and [ixe2, isel2[
                iex1 = j3 - nexcl2
                iex2 = j3 + nexcl2 + 1
                if nexcl == 0:
                    iex2 -= 1
                isel1 = max((j3 - nincl2), 0)
                isel2 = min((j3 + nincl2 + 1), nphibar)
                if isel1 == 0:
                    isel2 = min(nincl, nphibar)
                if isel2 == nphibar:
                    isel1 = max(nphibar - nincl, 0)

                # Loop over all Chi/Psi pixels
                for i2 in range(npsi):
                    for i1 in range(nchi):

                        # Apply Phibar sum correction
                        ipixu = j3*npix + i2*nchi + i1
                        if drb[ipixu] != 0.0:
                            sum_num = 0.0
                            if iex1 >= 1:
                                for i3 in range(isel1, iex1):
                                    ipix     = i3*npix + i2*nchi + i1
                                    sum_num += scrat[ipix]
                            if iex2 < nphibar:
                                for i3 in range(iex2, isel2):
                                    ipix     = i3*npix + i2*nchi + i1
                                    sum_num += scrat[ipix]
                            drb[ipixu] *= sum_num

            # Get statistics
            sum_dre = 0.0
            sum_drm = 0.0
            sum_drb = 0.0
            for i in range(dre.size()):
                sum_dre += dre[i]
                sum_drm += drm[i]
                sum_drb += drb[i]
            self._log_value(gammalib.NORMAL, 'Events in DRE', sum_dre)
            self._log_value(gammalib.NORMAL, 'Events in DRM', sum_drm)
            self._log_value(gammalib.NORMAL, 'Events in DRB', sum_drb)

            # Save DRB
            drb.save(drbfile, self['clobber'].boolean())

            # Log creation
            self._log_value(gammalib.NORMAL, 'DRB file created', drbfile.url())

        # Return
        return drbname

    def _generate_drb_bgdlixE(self, drename, drgname, drm, debug=False):
        """
        Generate DRB from DRG and DRE file with BGDLIX method E

        The method generates a background model using the BGDLIX method E,
        taking into account a source model DRM so that the method can be
        used for an iterative SRCLIX approach.

        Parameters
        ----------
        drename : str
            DRE filename
        drgname : str
            DRG filename
        drm : `~gammalib.GCOMDri`
            Source DRM
        debug : bool
            Write some files for debugging (default: False)

        Returns
        -------
        drbname : str
            DRB filename
        """
        # Get suffix and outfolder
        suffix    = self['suffix'].string()
        outfolder = self['outfolder'].string()

        # Get task parameters
        navgr = self['navgr'].integer()
        nincl = self['nincl'].integer()
        nexcl = self['nexcl'].integer()

        # Set DRB filename
        drbname = '%s/%s' % (outfolder, os.path.basename(drename.replace('dre',
                             'drb-bgdlixE-na%d-ni%d-ne%d%s' %
                             (navgr, nincl, nexcl, suffix))))
        drbfile = gammalib.GFilename(drbname)

        # Write header
        self._log_header3(gammalib.NORMAL, 'Compute DRB using BGDLIX method E')
        self._log_value(gammalib.NORMAL, 'Bins used for averaging', navgr)
        self._log_value(gammalib.NORMAL, 'Phibar layers to include', nincl)
        self._log_value(gammalib.NORMAL, 'Phibar layers to exclude', nexcl)

        # If DRB file exists then do nothing
        if drbfile.exists():
            self._log_value(gammalib.NORMAL, 'DRB file exists', drbfile.url())

        # ... otherwise compute and save it
        else:

            # Set DRE and DRG filenames
            drefile = gammalib.GFilename(drename)
            drgfile = gammalib.GFilename(drgname)

            # Load DRE and DRG
            dre = gammalib.GCOMDri(drefile)
            drg = gammalib.GCOMDri(drgfile)

            # Get dataspace dimensions
            nchi    = dre.nchi()
            npsi    = dre.npsi()
            nphibar = dre.nphibar()
            npix    = nchi * npsi

            # Multiply DRG pixels by solid angle
            for i in range(npix):
                omega = drg.map().solidangle(i)
                for i2 in range(nphibar):
                    index       = i + i2 * npix
                    drg[index] *= omega

            # Precompute half running average lengths
            navgr2 = int(navgr/2)
            nexcl2 = int(nexcl/2)
            nincl2 = int(nincl/2)

            # Initialise SCRAT and DRB
            scrat = drg.copy()
            drb   = drg.copy()

            # Initialise Chi2 and nincl2 map for bookkeeping
            chi2_map   = drg.map().copy()
            nincl2_map = drg.map().copy()
            chi2_map.stack_maps()
            nincl2_map.stack_maps()

            # Do Phibar normalization of DRG and store result in SCRAT.
            for i3 in range(nphibar):
                sum_dre = 0.0
                sum_drg = 0.0
                sum_drm = 0.0
                for i2 in range(npsi):
                    for i1 in range(nchi):
                        ipix     = i3*npix + i2*nchi + i1
                        sum_dre += dre[ipix]
                        sum_drg += drg[ipix]
                        sum_drm += drm[ipix]
                if sum_drg > 0.0:
                    scale = (sum_dre-sum_drm) / sum_drg
                    for i2 in range(npsi):
                        for i1 in range(nchi):
                            ipix        = i3*npix + i2*nchi + i1
                            scrat[ipix] = drg[ipix] * scale

            # Debug: save B^0_L
            if debug:
                srcatname = drbname.replace('drb', 'scrat1')
                srcatfile = gammalib.GFilename(srcatname)
                scrat.save(srcatfile, self['clobber'].boolean())

            # Initialise array of
            nincl2_hist = [0.0 for i in range(nincl2+1)]

            # Make background model
            for i1 in range(nchi):
                for i2 in range(npsi):

                    # Initialise nincl2 (add one since we reduce at beginning)
                    nincl2 = int(nincl/2) + 1

                    # Loop while nincl2 is larger than 1
                    while nincl2 > 1:

                        # Reduce nincl2
                        nincl2 -= 1

                        # Initialise Chi-squared
                        chi2 = 0.0
                        ndof = 0.0

                        # Loop over Phibar
                        for i3 in range(nphibar):

                            # Determine index range for sum. There are two sums that
                            # go over [isel1, iex1[ and [ixe2, isel2[
                            iex1 = i3 - nexcl2
                            iex2 = i3 + nexcl2 + 1
                            if nexcl == 0:
                                iex2 -= 1
                            isel1 = max((i3 - nincl2), 0)
                            isel2 = min((i3 + nincl2 + 1), nphibar)
                            if isel1 == 0:
                                isel2 = min(nincl, nphibar)
                            if isel2 == nphibar:
                                isel1 = max(nphibar - nincl, 0)

                            # Initialise sums
                            sum_dre   = 0.0
                            sum_drm   = 0.0
                            sum_scrat = 0.0

                            # Compute average over small patch in Chi/Psi
                            for j1 in range(max(0, i1 - navgr2), min(nchi, i1 + navgr2 + 1)):
                                for j2 in range(max(0, i2 - navgr2), min(npsi, i2 + navgr2 + 1)):

                                    # Take sum over [isel1, iex1[
                                    if iex1 >= 1:
                                        for j3 in range(isel1, iex1):
                                            jpix       = j1 + j2*nchi + j3*npix
                                            sum_dre   += dre[jpix]
                                            sum_drm   += drm[jpix]
                                            sum_scrat += scrat[jpix]

                                    # Take sum over [iex2, isel2[
                                    if iex2 < nphibar:
                                        for j3 in range(iex2, isel2):
                                            jpix       = j1 + j2*nchi + j3*npix
                                            sum_dre   += dre[jpix]
                                            sum_drm   += drm[jpix]
                                            sum_scrat += scrat[jpix]

                            # Re-normalise SCRAT to derive DRB
                            ipix = i1 + i2*nchi + i3*npix
                            if sum_scrat > 0.0:
                                drb[ipix] = scrat[ipix] * (sum_dre-sum_drm)/sum_scrat
                            else:
                                drb[ipix] = 0.0

                            # Compute Chi2
                            if drb[ipix] > 0.0:
                                arg   = dre[ipix] - drm[ipix] - drb[ipix]
                                chi2 += arg * arg / drb[ipix]
                                ndof += 1.0

                        # Normalise Chi2
                        if ndof > 0.0:
                            chi2 /= ndof

                        # Update maps
                        imap             = i1 + i2*nchi
                        chi2_map[imap]   = chi2
                        nincl2_map[imap] = nincl2

                        # Break if Chi2 is below threshold
                        if chi2 < 1.1:
                            break

                    # Log normalised Chi2
                    #msg   = '(Chi,Psi)=(%3d,%3d)' % (i1, i2)
                    #value = '%.2f (%d)' % (chi2, nincl2)
                    #self._log_value(gammalib.NORMAL, msg, value)

                    # Update nincl2 histogram
                    nincl2_hist[nincl2] += 1

            # Log nincl2 histogram
            for i in range(nincl/2+1):
                msg   = '(Chi,Psi) with nincl2 = %d' % (i)
                self._log_value(gammalib.NORMAL, msg, nincl2_hist[i])

            # Get statistics
            sum_dre = 0.0
            sum_drm = 0.0
            sum_drb = 0.0
            for i in range(dre.size()):
                sum_dre += dre[i]
                sum_drm += drm[i]
                sum_drb += drb[i]
            self._log_value(gammalib.NORMAL, 'Events in DRE', sum_dre)
            self._log_value(gammalib.NORMAL, 'Events in DRM', sum_drm)
            self._log_value(gammalib.NORMAL, 'Events in DRB', sum_drb)

            # Save DRB
            fits       = gammalib.GFits()
            drb_hdu    = drb.write(fits)
            chi2_hdu   = chi2_map.write(fits)
            nincl2_hdu = nincl2_map.write(fits)
            chi2_hdu.extname('Chi-squared')
            nincl2_hdu.extname('nincl2')
            fits.saveto(drbfile, self['clobber'].boolean())

            # Log creation
            self._log_value(gammalib.NORMAL, 'DRB file created', drbfile.url())

        # Return
        return drbname


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

        # Create database directory
        try:
            os.makedirs(gammalib.expand_env(self['outfolder'].string()))
        except OSError:
            pass

        # Log header
        self._log_header1(gammalib.NORMAL, 'Input observations')

        # Log input observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

        # Write header
        self._log_header1(gammalib.NORMAL, 'Extract models')

        # Extact extract source components that need to be subtracted from
        # input model
        models = gammalib.GModels()
        for model in self.obs().models():
            if model.type() == 'DRBFitting':
                self._log_value(gammalib.NORMAL, 'Exclude model', model.name())
            else:
                self._log_value(gammalib.NORMAL, 'Include model', model.name())
                models.append(model)
        self._log_value(gammalib.NORMAL, 'Number of source models', models.size())

        # Write header
        self._log_header1(gammalib.NORMAL, 'Create background model')

        # Get background method
        method = self['bkgmethod'].string()

        # Write background method
        self._log_value(gammalib.NORMAL, 'Background method', method)

        # Loop over all input observations
        for obs in self.obs():

            # Write header
            self._log_header2(gammalib.NORMAL, self._get_obs_header(obs))

            # Get DRE and DRG names
            drename = obs.drename().url()
            drgname = obs.drgname().url()

            # Compute DRM
            drm = obs.drm(models)

            # Generate background
            if method == 'PHINOR':
                drbname = self._generate_drb_phinor(drename, drgname, drm)
            elif method == 'BGDLIXA':
                drbname = self._generate_drb_bgdlixA(drename, drgname, drm)
            elif method == 'BGDLIXE':
                drbname = self._generate_drb_bgdlixE(drename, drgname, drm)

            # Set DRB name
            obs.drbname(drbname)

        # Log output observations
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
            # Log filename
            self._log_value(gammalib.NORMAL, 'Obs. definition XML file',
                                             outobs.url())

            # Save observations
            self.obs().save(outobs)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comobsback(sys.argv)

    # Execute application
    app.execute()
