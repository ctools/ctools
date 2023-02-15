#! /usr/bin/env python
# ==========================================================================
# Add COMPTEL observations
#
# Copyright (C) 2021-2023 Juergen Knoedlseder
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


# =============== #
# comobsadd class #
# =============== #
class comobsadd(ctools.csobservation):
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
        self._xml    = gammalib.GXml()
        self._dres   = []
        self._drbs   = []
        self._drws   = []
        self._drg    = []
        self._drx    = []
        self._nspuse = 0

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
        self['prefix'].string()
        self['coordsys'].string()
        self['proj'].string()
        if self['coordsys'].string() == 'CEL':
            self['ra'].real()
            self['dec'].real()
        else:
            self['glon'].real()
            self['glat'].real()
        self['nchi'].integer()
        self['npsi'].integer()
        self['dchi'].real()
        self['dpsi'].real()

        # Query ahead output model filename
        if self._read_ahead():
            self['outfolder'].string()
            self['outobs'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _create_obs(self, obs):
        """
        Create observation

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Observations

        Returns
        -------
        ebins : list of dict
            List with input DRIs
        """
        # Write header
        self._log_header1(gammalib.NORMAL, 'Create combined observation')

        # Create XML file
        self._xml = gammalib.GXml()

        # Initialise list of energy bins
        ebins = []

        # Continue only if there are observations
        if obs.size() > 0:

            # Append observation list to XML file
            xml_list = self._xml.append('observation_list title="observation list"')

            # Get dataspace dimensions
            dchi    = self['dchi'].real()
            dpsi    = self['dpsi'].real()
            dphibar = obs[0].events().dre().phibin()
            nchi    = self['nchi'].integer()
            npsi    = self['npsi'].integer()
            nphibar = obs[0].events().dre().nphibar()

            # Get centre of DRI
            coordsys = self['coordsys'].string()
            if self['coordsys'].string() == 'CEL':
                chi0 = self['ra'].real()
                psi0 = self['dec'].real()
            else:
                chi0 = self['glon'].real()
                psi0 = self['glat'].real()

            # Collect energy bins
            for o in obs:
                emin  = o.events().dre().ebounds().emin()
                emax  = o.events().dre().ebounds().emax()
                found = False
                for ebin in ebins:
                    if ebin['emin'] == emin and ebin['emax'] == emax:
                        ebin['dres'].append(o.drename())
                        ebin['drbs'].append(o.drbname())
                        ebin['drws'].append(o.drwname())
                        ebin['drgs'].append(o.drgname())
                        ebin['drxs'].append(o.drxname())
                        ebin['gti'].extend(o.events().dre().gti())
                        ebin['nspinp'] += o.events().dre().num_superpackets()
                        ebin['nspuse'] += o.events().dre().num_used_superpackets()
                        ebin['nspskp'] += o.events().dre().num_skipped_superpackets()
                        found = True
                        break
                if not found:
                    ebins.append({'emin'     : emin,
                                  'emax'     : emax,
                                  'dres'     : [o.drename()],
                                  'drbs'     : [o.drbname()],
                                  'drws'     : [o.drwname()],
                                  'drgs'     : [o.drgname()],
                                  'drxs'     : [o.drxname()],
                                  'iaq'      : o.response().rspname(),
                                  'gti'      : o.events().dre().gti(),
                                  'tofcor'   : o.events().dre().tof_correction(),
                                  'phasecor' : o.events().dre().phase_correction(),
                                  'nspinp'   : o.events().dre().num_superpackets(),
                                  'nspuse'   : o.events().dre().num_used_superpackets(),
                                  'nspskp'   : o.events().dre().num_skipped_superpackets()})

            # Log DRI definition
            self._log_value(gammalib.NORMAL, 'Coordinate system', coordsys)
            self._log_value(gammalib.NORMAL, 'Chi0 (deg)', chi0)
            self._log_value(gammalib.NORMAL, 'Psi0 (deg)', psi0)
            self._log_value(gammalib.NORMAL, 'Chi binsize (deg)', dchi)
            self._log_value(gammalib.NORMAL, 'Psi binsize (deg)', dpsi)
            self._log_value(gammalib.NORMAL, 'Phibar binsize (deg)', dphibar)
            self._log_value(gammalib.NORMAL, 'Number of Chi bins', nchi)
            self._log_value(gammalib.NORMAL, 'Number of Psi bins', npsi)
            self._log_value(gammalib.NORMAL, 'Number of Phibar bins', nphibar)
            self._log_value(gammalib.NORMAL, 'Number of energy bins', len(ebins))
            for i, ebin in enumerate(ebins):
                self._log_value(gammalib.NORMAL, 'Energy bin %d' % i,
                                '%s - %s' % (str(ebin['emin']), str(ebin['emax'])))
                self._log_value(gammalib.NORMAL, ' ToF correction', ebin['tofcor'])
                self._log_value(gammalib.NORMAL, ' Phase correction', ebin['phasecor'])
                self._log_value(gammalib.NORMAL, ' Response', ebin['iaq'])

            # Create sky map for DRE, DRB, DRW and DRG
            cube = gammalib.GSkyMap(self['proj'].string(), coordsys,
                                    chi0, psi0, -dchi, dpsi, nchi, npsi)

            # Create sky map for DRX
            expo = gammalib.GSkyMap('CAR', 'GAL', 0.0, 0.0, -1.0, 1.0, 360, 180)

            # Allocate DRE, DRB, DRW, DRG and DRX
            dre       = gammalib.GCOMDri(cube, 0.0, dphibar, nphibar)
            drb       = gammalib.GCOMDri(cube, 0.0, dphibar, nphibar)
            drw       = gammalib.GCOMDri(cube, 0.0, dphibar, nphibar)
            self._drg = gammalib.GCOMDri(cube, 0.0, dphibar, nphibar)
            self._drx = gammalib.GCOMDri(expo)

            # Set GTIs of DRE, DRB, DRW, DRG and DRX
            dre.gti(ebin['gti'])
            drb.gti(ebin['gti'])
            drw.gti(ebin['gti'])
            self._drg.gti(ebin['gti'])
            self._drx.gti(ebin['gti'])

            # Set DRG and DRX filenames
            prefix = self['prefix'].string()
            self._drg.name('%s/%s_drg.fits' % (self['outfolder'].string(), prefix))
            self._drx.name('%s/%s_drx.fits' % (self['outfolder'].string(), prefix))

            # Set superpackets usage
            if len(ebins) > 0:
                self._drg.num_superpackets(ebins[0]['nspinp'])
                self._drg.num_used_superpackets(ebins[0]['nspuse'])
                self._drg.num_skipped_superpackets(ebins[0]['nspskp'])
                self._drx.num_superpackets(ebins[0]['nspinp'])
                self._drx.num_used_superpackets(ebins[0]['nspuse'])
                self._drx.num_skipped_superpackets(ebins[0]['nspskp'])

            # Create DRE, DRB and DRW for each energy bin and append an XML entry
            for ebin in ebins:

                # Set observation name and id
                name = 'Combined'
                id   = '%6.6d-%6.6dkeV' % (ebin['emin'].keV(), ebin['emax'].keV())

                # Set energy boundaries of DRE and DRB
                ebounds = gammalib.GEbounds()
                ebounds.append(ebin['emin'], ebin['emax'])
                dre.ebounds(ebounds)
                drb.ebounds(ebounds)
                drw.ebounds(ebounds)

                # Set ToF correction
                dre.tof_correction(ebin['tofcor'])

                # Set phase correction
                dre.phase_correction(ebin['phasecor'])

                # Set superpackets usage
                dre.num_superpackets(ebin['nspinp'])
                dre.num_used_superpackets(ebin['nspuse'])
                dre.num_skipped_superpackets(ebin['nspskp'])
                drb.num_superpackets(ebin['nspinp'])
                drb.num_used_superpackets(ebin['nspuse'])
                drb.num_skipped_superpackets(ebin['nspskp'])
                drw.num_superpackets(ebin['nspinp'])
                drw.num_used_superpackets(ebin['nspuse'])
                drw.num_skipped_superpackets(ebin['nspskp'])

                # Set DRE, DRB and DRW filenames
                dre.name('%s/%s_%s_dre.fits' % (self['outfolder'].string(), prefix, id))
                drb.name('%s/%s_%s_drb.fits' % (self['outfolder'].string(), prefix, id))
                drw.name('%s/%s_%s_drw.fits' % (self['outfolder'].string(), prefix, id))

                # Put DRE, DRB and DRW in lists
                self._dres.append(dre.copy())
                self._drbs.append(drb.copy())
                self._drws.append(drw.copy())

                # Append observation
                xml_obs = xml_list.append('observation name="%s" id="%s" instrument="COM"' % \
                                          (name, id))

                # Append DRE, DRB, DRW, DRG and DRX files
                xml_obs.append('parameter name="DRE" file="%s"' % (dre.name()))
                xml_obs.append('parameter name="DRB" file="%s"' % (drb.name()))
                xml_obs.append('parameter name="DRW" file="%s"' % (drw.name()))
                xml_obs.append('parameter name="DRG" file="%s"' % (self._drg.name()))
                xml_obs.append('parameter name="DRX" file="%s"' % (self._drx.name()))

                # Append IAQ file
                xml_obs.append('parameter name="IAQ" value="%s"' % (ebin['iaq']))

        # Return energy bins
        return ebins

    def _combine_obs(self, ebins):
        """
        Combine observations

        Parameters
        ----------
        ebins : list of dict
            List with input DRIs
        """
        # Write header
        self._log_header1(gammalib.NORMAL, 'Combine observations')

        # Combine DREs
        self._combine_dres(ebins)

        # Combine DRBs
        self._combine_drbs(ebins)

        # Combine DRWs
        self._combine_drws(ebins)

        # Combine DRGs
        self._combine_drgs(ebins)

        # Combine DRXs
        self._combine_drxs(ebins)

        # Return
        return

    def _combine_dres(self, ebins):
        """
        Combine DREs

        The DREs are combined by redistributing the individual events within
        their Chi/Psi pixels using a random number generator. In that way no
        precise match between the original bins and the bins of the result
        DRE is needed.

        Parameters
        ----------
        ebins : list of dict
            List with input DRIs
        """
        # Write header
        self._log_header2(gammalib.NORMAL, 'Combine DREs')

        # Initialise random number generator
        ran = gammalib.GRan()

        # Loop over all energy bins
        for iebin, ebin in enumerate(ebins):

            # Write header
            self._log_header3(gammalib.NORMAL, 'Energy bin %s - %s' %
                              (str(ebin['emin']), str(ebin['emax'])))

            # Set combined data space map
            map = self._dres[iebin].map()

            # Loop over DREs
            for drename in ebin['dres']:

                # Log file name
                self._log_value(gammalib.NORMAL, 'Filename', drename.url())

                # Load DRE
                dre = gammalib.GCOMDri(drename)

                # Get dimensions of data space
                npix    = dre.nchi() * dre.npsi()
                nphibar = dre.nphibar()

                # Initialise number of distributed events
                num_events  = 0.0
                num_sampled = 0
                num_added   = 0

                # Loop over data space bins
                for iphibar in range(nphibar):
                    for ipix in range(npix):

                        # Get number of events
                        events      = dre.map()(ipix, iphibar)
                        num_events += events

                        # Get pixel in original DRE
                        pixel_src = dre.map().inx2pix(ipix)

                        # Loop over all events
                        for i in range(int(events+0.5)):

                            # Randomise sky pixel
                            x     = pixel_src.x() + (ran.uniform() - 0.5)
                            y     = pixel_src.y() + (ran.uniform() - 0.5)
                            pixel = gammalib.GSkyPixel(x, y)

                            # Get sky direction for randomised sky pixel
                            dir = dre.map().pix2dir(pixel)

                            # Put event in combined data space
                            if map.contains(dir):
                                idst                = map.dir2inx(dir)
                                map[idst, iphibar] += 1.0
                                num_added          += 1

                            # Update number of sampled events
                            num_sampled += 1

                # Log number of events
                self._log_value(gammalib.NORMAL, 'Events in DRE', num_events)
                self._log_value(gammalib.NORMAL, 'Re-distributed events', num_sampled)
                self._log_value(gammalib.NORMAL, 'Events inside added DRE', num_added)
                self._log_value(gammalib.NORMAL, 'Events outside added DRE', num_sampled-num_added)

        # Return
        return

    def _combine_drbs(self, ebins):
        """
        Combine DRBs

        The DRBs are combined by dividing the original pixels by their solid
        angles and multipliying with the solid angle of the combined DRB. In
        that way the solid angle changes due to different projections are
        correctly taken into account.

        Parameters
        ----------
        ebins : list of dict
            List with input DRIs
        """
        # Write header
        self._log_header2(gammalib.NORMAL, 'Combine DRBs')

        # Loop over all energy bins
        for iebin, ebin in enumerate(ebins):

            # Write header
            self._log_header3(gammalib.NORMAL, 'Energy bin %s - %s' %
                              (str(ebin['emin']), str(ebin['emax'])))

            # Get dimensions of combined data space
            npix    = self._drbs[iebin].nchi() * self._drbs[iebin].npsi()
            nphibar = self._drbs[iebin].nphibar()

            # Set combined data space map
            map = self._drbs[iebin].map()

            # Loop over DRBs
            for drbname in ebin['drbs']:

                # Log DRB file name
                self._log_value(gammalib.NORMAL, 'Filename', drbname.url())

                # Load DRB
                drb = gammalib.GCOMDri(drbname)

                # Divide pixels by solid angle so that the DRB returns the
                # number of background counts per steradian
                nmaps   = drb.map().nmaps()
                npixels = drb.map().npix()
                for ipixel in range(npixels):
                    omega = drb.map().solidangle(ipixel)
                    for imap in range(nmaps):
                        drb.map()[ipixel,imap] /= omega

                # Loop over combined data space Chi/Psi pixels
                for ipix in range(npix):

                    # Get Chi/Psi direction of pixel
                    dir = map.inx2dir(ipix)

                    # If Chi/Psi is contained in DRB then extract value.
                    # Put the check in a try-except block since there
                    # may be coordinate transformation issues at this
                    # point
                    try:
                        if drb.map().contains(dir):

                            # Get solid angle of target map to convert the
                            # DRB values into absolute counts
                            omega = map.solidangle(ipix)

                            # Loop over combined data space Phibar layers
                            for iphibar in range(nphibar):
                                map[ipix, iphibar] += drb.map()(dir, iphibar) * omega

                    except RuntimeError:
                        pass

        # Return
        return

    def _combine_drws(self, ebins):
        """
        Combine DRWs

        The DRWs are combined by dividing the original pixels by their solid
        angles and multipliying with the solid angle of the combined DRW. In
        that way the solid angle changes due to different projections are
        correctly taken into account.

        Parameters
        ----------
        ebins : list of dict
            List with input DRIs
        """
        # Write header
        self._log_header2(gammalib.NORMAL, 'Combine DRWs')

        # Loop over all energy bins
        for iebin, ebin in enumerate(ebins):

            # Write header
            self._log_header3(gammalib.NORMAL, 'Energy bin %s - %s' %
                              (str(ebin['emin']), str(ebin['emax'])))

            # Get dimensions of combined data space
            npix    = self._drws[iebin].nchi() * self._drws[iebin].npsi()
            nphibar = self._drws[iebin].nphibar()

            # Set combined data space map
            map = self._drws[iebin].map()

            # Loop over DRWs
            for drwname in ebin['drws']:

                # Log DRW file name
                self._log_value(gammalib.NORMAL, 'Filename', drwname.url())

                # Load DRW
                drw = gammalib.GCOMDri(drwname)

                # Divide pixels by solid angle so that the DRW returns the
                # number of background counts per steradian
                nmaps   = drw.map().nmaps()
                npixels = drw.map().npix()
                for ipixel in range(npixels):
                    omega = drw.map().solidangle(ipixel)
                    for imap in range(nmaps):
                        drw.map()[ipixel,imap] /= omega

                # Loop over combined data space Chi/Psi pixels
                for ipix in range(npix):

                    # Get Chi/Psi direction of pixel
                    dir = map.inx2dir(ipix)

                    # If Chi/Psi is contained in DRW then extract value.
                    # Put the check in a try-except block since there
                    # may be coordinate transformation issues at this
                    # point
                    try:
                        if drw.map().contains(dir):

                            # Get solid angle of target map to convert the
                            # DRB values into absolute counts
                            omega = map.solidangle(ipix)

                            # Loop over combined data space Phibar layers
                            for iphibar in range(nphibar):
                                map[ipix, iphibar] += drw.map()(dir, iphibar) * omega

                    except RuntimeError:
                        pass

        # Return
        return

    def _combine_drgs(self, ebins):
        """
        Combine DRGs

        Parameters
        ----------
        ebins : list of dict
            List with input DRIs
        """
        # Write header
        self._log_header2(gammalib.NORMAL, 'Combine DRGs')

        # Continue only if there are energy bins
        if len(ebins) > 0:

            # Get dimensions of combined data space
            npix    = self._drg.nchi() * self._drg.npsi()
            nphibar = self._drg.nphibar()

            # Set combined data space map
            map = self._drg.map()

            # Loop over all DRGs
            for drgname in ebins[0]['drgs']:

                # Log DRG file name
                self._log_header3(gammalib.NORMAL, drgname.url())

                # Load DRG
                drg = gammalib.GCOMDri(drgname)

                # Compute DRG weight
                weight = float(drg.num_used_superpackets()) / float(ebins[0]['nspuse'])
                self._log_value(gammalib.NORMAL, 'Weight', weight)

                # Loop over combined data space Chi/Psi pixels
                for ipix in range(npix):

                    # Get Chi/Psi direction of pixel
                    dir = map.inx2dir(ipix)

                    # If Chi/Psi is contained in DRB then extract value.
                    # Put the check in a try-except block since there
                    # may be coordinate transformation issues at this
                    # point
                    try:
                        if drg.map().contains(dir):

                            # Loop over combined data space Phibar layers
                            for iphibar in range(nphibar):
                                map[ipix, iphibar] += weight * drg.map()(dir, iphibar)
                    except RuntimeError:
                        pass

        # Return
        return

    def _combine_drxs(self, ebins):
        """
        Combine DRXs

        Parameters
        ----------
        ebins : list of dict
            List with input DRIs
        """
        # Write header
        self._log_header2(gammalib.NORMAL, 'Combine DRXs')

        # Continue only if there are energy bins
        if len(ebins) > 0:

            # Initialise DRX value
            drx_value = 0.0

            # Loop over all DRXs
            for drxname in ebins[0]['drxs']:

                # Log DRX file name
                self._log_header3(gammalib.NORMAL, drxname.url())

                # Load DRX
                drx = gammalib.GCOMDri(drxname)

                # Compute DRX maximum
                drx_max = 0.0
                for ipix in range(drx.size()):
                    value = drx[ipix]
                    if value > drx_max:
                        drx_max = value

                # Add maximum to DRX value
                drx_value += drx_max

                # Write maximum into logger
                self._log_value(gammalib.NORMAL, 'Maximum', '%.1f cm^2 s' % (drx_max))

            # Write maximum into logger
            self._log_value(gammalib.NORMAL, 'Total DRX value', '%.1f cm^2 s' % (drx_value))

            # Set combined DRX
            for ipix in range(self._drx.size()):
                self._drx[ipix] = drx_value

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
        self._log_header1(gammalib.NORMAL, 'Input observations')

        # Log input observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

        # Create combined observation
        ebins = self._create_obs(self.obs())

        # Combine observations
        self._combine_obs(ebins)

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

            # Create outfolder directory
            try:
                os.makedirs(gammalib.expand_env(self['outfolder'].string()))
            except OSError:
                pass

            # Save observation definition file
            self._xml.save(outobs)
            self._log_value(gammalib.NORMAL, 'Obs. definition XML file', outobs.url())

            # Save DREs
            for dre in self._dres:
                dre.save(dre.name(), self['clobber'].boolean())
                self._log_value(gammalib.NORMAL, 'DRE file', dre.name())

            # Save DRBs
            for drb in self._drbs:
                drb.save(drb.name(), self['clobber'].boolean())
                self._log_value(gammalib.NORMAL, 'DRB file', drb.name())

            # Save DRWs
            for drw in self._drws:
                drw.save(drw.name(), self['clobber'].boolean())
                self._log_value(gammalib.NORMAL, 'DRW file', drw.name())

            # Save DRG
            self._drg.save(self._drg.name(), self['clobber'].boolean())
            self._log_value(gammalib.NORMAL, 'DRG file', self._drg.name())

            # Save DRX
            self._drx.save(self._drx.name(), self['clobber'].boolean())
            self._log_value(gammalib.NORMAL, 'DRX file', self._drx.name())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = comobsadd(sys.argv)

    # Execute application
    app.execute()
