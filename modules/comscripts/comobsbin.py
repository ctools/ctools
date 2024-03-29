#! /usr/bin/env python
# ==========================================================================
# Bin COMPTEL observations
#
# Copyright (C) 2019-2023 Juergen Knoedlseder
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
# comobsbin class #
# =============== #
class comobsbin(ctools.csobservation):
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
        self._phases           = gammalib.GPhases()
        self._select           = gammalib.GCOMSelection()
        self._drx_suffix       = ''
        self._drg_suffix       = ''
        self._drw_suffix       = ''
        self._dre_suffix       = ''
        self._global_datastore = '$COMDATA/datastore'

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

        # Query parameters
        self['response'].string()
        self._get_ebounds()
        self._get_phases()
        if not self._phases.is_empty():
            if gammalib.toupper(self['psrname'].string()) != 'NONE':
                self['ephemerides'].filename()
            else:
                self['phase0'].time()
                self['period'].real()
        self['coordsys'].string()
        self['proj'].string()
        if not self['usepnt'].boolean():
            self['chi0'].real()
            self['psi0'].real()
        self['outfolder'].string()
        self['dchi'].real()
        self['dpsi'].real()
        self['dphibar'].real()
        self['nchi'].integer()
        self['npsi'].integer()
        self['nphibar'].integer()
        self['psdmin'].integer()
        self['psdmax'].integer()
        self['zetamin'].real()
        self['fpmtflag'].integer()
        self['drwmethod'].string()
        self['timebin'].real()

        # Get D1 and D2 module usage strings
        d1use = self['d1use'].string()
        d2use = self['d2use'].string()

        # Check D1 and D2 module usage strings
        if len(d1use) != 7:
            msg = 'Incorrect length %d of D1 usage string. String needs to have 7 digits.' % \
                  len(d1use)
            raise RuntimeError(msg)
        if len(d2use) != 14:
            msg = 'Incorrect length %d of D2 usage string.  String needs to have 14 digits.' % \
                  len(d2use)
            raise RuntimeError(msg)

        # Query ahead output model filename
        if self._read_ahead():
            self['outobs'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _get_ebounds(self):
        """
        Get energy boundaries
        """
        # Get type of energy boundaries
        ebinalg = self['ebinalg'].string()

        # Handle energy boundaries
        if ebinalg == 'FILE':
            ebounds = gammalib.GEbounds(self['ebinfile'].filename())
        else:
            emin     = gammalib.GEnergy(self['emin'].real(), 'MeV')
            emax     = gammalib.GEnergy(self['emax'].real(), 'MeV')
            enumbins = self['enumbins'].integer()
            ebounds  = gammalib.GEbounds(enumbins, emin, emax, ebinalg)

        # Return energy boundaries
        return ebounds

    def _get_phases(self):
        """
        Get phases
        """
        # Initialise phases
        self._phases.clear()

        # Continue only if phase string is valid
        if gammalib.toupper(self['phase'].string()) != 'NONE':

            # Get phases and split them
            phases = self['phase'].string().split(';')

            # Loop over phases
            for phase in phases:

                # Get phase boundaries
                phase_bounds = phase.split('-')
                if len(phase_bounds) != 2:
                    msg = 'Invalid phase selection string "'+phase+'" encountered. '+\
                          'The phase selection string requires a minimum and '+\
                          'maximum value separated by a colon, e.g. "0.3-0.5".'
                    raise RuntimeError(msg)
                elif len(phase_bounds[0]) == 0 or len(phase_bounds[1]) == 0:
                    msg = 'Invalid phase selection string "'+phase+'" encountered. '+\
                          'The phase selection string requires a minimum and '+\
                          'maximum value separated by a colon, e.g. "0.3-0.5".'
                    raise RuntimeError(msg)
                pmin = float(phase_bounds[0])
                pmax = float(phase_bounds[1])

                # Check phase boundaries
                if pmin < 0.0 or pmin > 1.0:
                    msg = 'Phase minimum %f outside the valid range [0,1]. Please '+\
                          'specify phase interval boundaries comprised within 0 '+\
                          'and 1.' % pmin
                    raise RuntimeError(msg)
                if pmax < 0.0 or pmax > 1.0:
                    msg = 'Phase maximum %f outside the valid range [0,1]. Please '+\
                          'specify phase interval boundaries comprised within 0 '+\
                          'and 1.' % pmax
                    raise RuntimeError(msg)

                # Append phase interval
                self._phases.append(pmin, pmax)

        # Return
        return

    def _set_selection_set(self):
        """
        Set selection set
        """
        # Write header
        self._log_header1(gammalib.NORMAL, 'Set selection set')

        # Clear selection set
        self._select.clear()

        # If phases are set then handle phase selection
        if not self._phases.is_empty():

            # If pulsar string is valid then handle pulsar
            psrname = self['psrname'].string()
            if gammalib.toupper(psrname) != 'NONE':

                # Log pulsar phase selection
                self._log_string(gammalib.NORMAL, 'Phase selection for pulsar "%s"' %
                                 psrname)

                # Read ephemerides
                pulsar = gammalib.GPulsar(self['ephemerides'].filename(), psrname)

                # Set pulsar ephemerides
                self._select.pulsar(pulsar)

                # Set pulsar phases
                self._select.pulsar_phases(self._phases)

                # Set phase suffix
                suffix = '_pulphase'

                # Set DRE, DRG and DRX suffix
                phases           = self['phase'].string()
                self._dre_suffix = suffix+phases.replace(';', '_')
                self._drg_suffix = suffix
                self._drw_suffix = suffix
                self._drx_suffix = suffix

            # ... otherwise handle orbital phase selection
            else:

                # Raise error if period is not positive
                period = self['period'].real()
                if period <= 0.0:
                    msg = 'Period %f is not positive. Please specify a positive period.'
                    raise RuntimeError(msg)

                # Log orbital phase selection
                self._log_string(gammalib.NORMAL, 'Orbital selection for period %f days' %
                                 period)

                # Set orbital period
                self._select.orbital_period(period, self['phase0'].time())

                # Set orbital phases
                self._select.orbital_phases(self._phases)

                # Set phase suffix
                suffix = '_orbphase'

                # Set DRE, DRG and DRX suffix
                phases           = self['phase'].string()
                self._dre_suffix = suffix+phases.replace(';', '_')
                self._drg_suffix = suffix+phases.replace(';', '_')
                self._drw_suffix = suffix+phases.replace(';', '_')
                self._drx_suffix = suffix+phases.replace(';', '_')

        # If PSD interval differs from standard interval then set the interval
        # and append flag to suffix
        if self['psdmin'].integer() != 0 or self['psdmax'].integer() != 110:

            # Set PSD interval
            self._select.psd_min(self['psdmin'].integer())
            self._select.psd_max(self['psdmax'].integer())

            # Set DRE suffix
            self._dre_suffix += '_psd%d-%d' % (self['psdmin'].integer(),
                                               self['psdmax'].integer())

        # If zeta angle differs from standard value then append zeta angle to DRE
        # DRG and DRW suffix
        if self['zetamin'].real() != 5.0:
            self._dre_suffix += '_zeta%.1f' % (self['zetamin'].real())
            self._drg_suffix += '_zeta%.1f' % (self['zetamin'].real())
            self._drw_suffix += '_zeta%.1f' % (self['zetamin'].real())

        # If handling of D2 modules with failed PMT flag differs from standard then
        # set the flag and append flag to suffix
        if self['fpmtflag'].integer() != 0:

            # Set handling of D2 modules with failed PMT flag
            self._select.fpmtflag(self['fpmtflag'].integer())

            # Set DRE, DRG and DRW suffix
            self._dre_suffix += '_fpmt%1d' % (self['fpmtflag'].integer())
            self._drg_suffix += '_fpmt%1d' % (self['fpmtflag'].integer())
            self._drw_suffix += '_fpmt%1d' % (self['fpmtflag'].integer())

        # If D1 module usage differs from standard value then set module usage
        # and append usage to DRE and DRG suffix
        d1use = self['d1use'].string()
        if d1use != '1111111':

            # Set usage flags
            for i in range(7):
                if d1use[i] == '1':
                    self._select.use_d1(i,True)
                else:
                    self._select.use_d1(i,False)

            # Set DRE, DRG and DRW suffix
            self._dre_suffix += '_%s' % (d1use)
            self._drg_suffix += '_%s' % (d1use)
            self._drw_suffix += '_%s' % (d1use)

        # If D2 module usage differs from standard value then set module usage
        # and append usage to DRE and DRG suffix
        d2use = self['d2use'].string()
        if d2use != '11111111111111':

            # Set usage flags
            for i in range(14):
                if d2use[i] == '1':
                    self._select.use_d2(i,True)
                else:
                    self._select.use_d2(i,False)

            # Set DRE, DRG and DRW suffix
            self._dre_suffix += '_%s' % (d2use)
            self._drg_suffix += '_%s' % (d2use)
            self._drw_suffix += '_%s' % (d2use)

        # Log selection set
        self._log_string(gammalib.NORMAL, str(self._select))

        # Return
        return

    def _get_pointing_direction(self, obs):
        """
        Get pointing direction from observation
        """
        # Get copy of pointing direction from ROI centre
        dir = obs.events().roi().centre().dir().copy()

        # Return pointing direction
        return dir

    def _generate_dri(self, obs, ebounds):
        """
        Generate DRI

        Parameters
        ----------
        obs : `~gammalib.GCOMObservation`
            COMPTEL observation
        ebounds : `~gammalib.GEbounds`
            Energy boundaries

        Returns
        -------
        drxname, drgname, drwnames, drenames : tuple of str
            Names of DRX, DRG, DRWs and DREs
        """
        # Get pointing direction
        pnt = self._get_pointing_direction(obs)

        # Get DRE and DRG centre definition, coordinate system and projection
        coordsys = self['coordsys'].string()
        proj     = self['proj'].string()
        if self['usepnt'].boolean():
            if coordsys == 'GAL':
                chi0 = pnt.l_deg()
                psi0 = pnt.b_deg()
            else:
                chi0 = pnt.ra_deg()
                psi0 = pnt.dec_deg()
        else:
            chi0 = self['chi0'].real()
            psi0 = self['psi0'].real()

        # Get DRE and DRG dimension
        dchi    = self['dchi'].real()
        dpsi    = self['dpsi'].real()
        dphibar = self['dphibar'].real()
        nchi    = self['nchi'].integer()
        npsi    = self['npsi'].integer()
        nphibar = self['nphibar'].integer()

        # Set DRI prefix in case that the DRI definition deviates from the standard
        dri_prefix = ''
        filler     = '_'
        if coordsys != 'GAL':
            dri_prefix += '_cel'
            filler      = '-'
        if proj != 'TAN':
            dri_prefix += '%s%s' % (filler, proj)
            filler      = '-'
        if nchi != 80 or npsi != 80 or nphibar != 25:
            dri_prefix += '%s%dx%dx%d' % (filler, nchi, npsi, nphibar)
            filler      = '-'
        if dchi != 1.0 or dpsi != 1.0 or dphibar != 2.0:
            dri_prefix += '%s%.1fx%.1fx%.1f' % (filler, dchi, dpsi, dphibar)
            filler      = '-'
        if not self['usepnt'].boolean():
            if psi0 < 0.0:
                dri_prefix += '%s%.2f-%.2f' % (filler, chi0,-psi0)
            else:
                dri_prefix += '%s%.2f+%.2f' % (filler, chi0,psi0)

        # Create sky map for DRE and DRG
        cube = gammalib.GSkyMap(proj, coordsys, chi0, psi0, dchi, dpsi, nchi, npsi)

        # Create sky map for DRX
        expo = gammalib.GSkyMap('CAR', 'GAL', 0.0, 0.0, 1.0, 1.0, 360, 180)

        # Allocate DRE, DRW, DRG and DRX
        dre = gammalib.GCOMDri(cube, 0.0, dphibar, nphibar)
        drw = gammalib.GCOMDri(cube, 0.0, dphibar, nphibar)
        drg = gammalib.GCOMDri(cube, 0.0, dphibar, nphibar)
        drx = gammalib.GCOMDri(expo)

        # Allocate container for DRWs and initialise filename list
        drws     = gammalib.GCOMDris()
        drwfiles = []

        # Set DRG and DRX filenames in output folder
        drxname = '%s/%s_drx%s.fits'   % (self['outfolder'].string(), obs.id(), self._drx_suffix)
        drgname = '%s/%s%s_drg%s.fits' % (self['outfolder'].string(), obs.id(), dri_prefix, self._drg_suffix)
        drxfile = gammalib.GFilename(drxname)
        drgfile = gammalib.GFilename(drgname)

        # Set DRG and DRX filenames in global data store
        drxname_global = '%s/%s_drx%s.fits'   % (self._global_datastore, obs.id(), self._drx_suffix)
        drgname_global = '%s/%s%s_drg%s.fits' % (self._global_datastore, obs.id(), dri_prefix, self._drg_suffix)
        drxfile_global = gammalib.GFilename(drxname_global)
        drgfile_global = gammalib.GFilename(drgname_global)

        # Compute DRX
        self._log_header3(gammalib.NORMAL, 'Compute DRX')
        if drxfile_global.exists():
            drxname = drxname_global
            self._log_value(gammalib.NORMAL, 'Global DRX file exists', drxfile_global.url())
        elif drxfile.exists():
            self._log_value(gammalib.NORMAL, 'Local DRX file exists', drxfile.url())
        else:
            drx.compute_drx(obs, self._select)
            drx.save(drxfile, True)
            self._log_value(gammalib.NORMAL, 'Local DRX file created', drxfile.url())
            self._log_string(gammalib.NORMAL, str(drx))

        # Compute DRG
        self._log_header3(gammalib.NORMAL, 'Compute DRG')
        if drgfile_global.exists():
            drgname = drgname_global
            self._log_value(gammalib.NORMAL, 'Global DRG file exists', drgfile_global.url())
        elif drgfile.exists():
            self._log_value(gammalib.NORMAL, 'Local DRG file exists', drgfile.url())
        else:
            drg.compute_drg(obs, self._select, self['zetamin'].real())
            drg.save(drgfile, True)
            self._log_value(gammalib.NORMAL, 'Local DRG file created', drgfile.url())
            self._log_string(gammalib.NORMAL, str(drg))

        # Initialse list of DRE names
        drenames = []

        # Generate one DRE for each energy boundary
        for i in range(ebounds.size()):

            # Set DRE filename in output folder
            drename = '%s/%s%s_dre%s_%6.6d-%6.6dkeV.fits' % \
                      (self['outfolder'].string(), obs.id(), dri_prefix, self._dre_suffix,
                       ebounds.emin(i).keV(), ebounds.emax(i).keV())
            drefile = gammalib.GFilename(drename)

            # Set DRE filename in global data store
            drename_global = '%s/%s%s_dre%s_%6.6d-%6.6dkeV.fits' % \
                             (self._global_datastore, obs.id(), dri_prefix, self._dre_suffix,
                              ebounds.emin(i).keV(), ebounds.emax(i).keV())
            drefile_global = gammalib.GFilename(drename_global)

            # Write header
            self._log_header3(gammalib.NORMAL, 'Compute DRE for %.3f - %.3f MeV' % \
                              (ebounds.emin(i).MeV(), ebounds.emax(i).MeV()))

            # If DRE file exists in global datastor then do nothing
            if drefile_global.exists():
                drename = drename_global
                self._log_value(gammalib.NORMAL, 'Global DRE file exists', drefile_global.url())

            # If DRE file exists in output folder then do nothing
            elif drefile.exists():
                self._log_value(gammalib.NORMAL, 'Local DRE file exists', drefile.url())

            # ... otherwise compute and save it
            else:

                # Set DRE energy range
                dre.ebounds(gammalib.GEbounds(ebounds.emin(i), ebounds.emax(i)))

                # Compute DRE
                dre.compute_dre(obs, self._select, self['zetamin'].real())

                # Save DRE
                dre.save(drefile, True)

                # Log creation
                self._log_value(gammalib.NORMAL, 'Local DRE file created', drefile.url())

                # Log DRE
                self._log_string(gammalib.NORMAL, str(dre))

            # Append DRE filename in output folder
            drenames.append(drename)

        # Initialse list of DRW names
        drwnames = []
        engindex = []

        # Get DRW method and suffix in lower case. Add a kludge that avoids appending
        # a suffix for the 'phibar' method with a timebin of 300 for backwards
        # compatibility with already existing DRW files.
        drwmethod = gammalib.tolower(self['drwmethod'].string())
        drwsuffix = '-%s' % (drwmethod)
        if drwmethod == 'phibar':
            drwsuffix += '%d' % (self['timebin'].real())
        if drwsuffix == '-phibar300':
            drwsuffix = ''

        # Generate one DRW for each energy boundary
        for i in range(ebounds.size()):

            # Set DRW filename in output folder
            drwname = '%s/%s%s_drw%s%s_%6.6d-%6.6dkeV.fits' % \
                       (self['outfolder'].string(), obs.id(), dri_prefix, drwsuffix,
                        self._drw_suffix, ebounds.emin(i).keV(), ebounds.emax(i).keV())
            drwfile = gammalib.GFilename(drwname)

            # Set DRW filename in global data store
            drwname_global = '%s/%s%s_drw%s%s_%6.6d-%6.6dkeV.fits' % \
                              (self._global_datastore, obs.id(), dri_prefix, drwsuffix,
                               self._drw_suffix, ebounds.emin(i).keV(), ebounds.emax(i).keV())
            drwfile_global = gammalib.GFilename(drwname_global)

            # Write header
            self._log_header3(gammalib.NORMAL, 'Check DRW for %.3f - %.3f MeV' % \
                                      (ebounds.emin(i).MeV(), ebounds.emax(i).MeV()))

            # If DRW file exists in global datastore then do nothing
            if drwfile_global.exists():
                drwname = drwname_global
                self._log_value(gammalib.NORMAL, 'Global DRW file exists', drwfile_global.url())

            # If DRW file exists in output folder then do nothing
            elif drwfile.exists():
                self._log_value(gammalib.NORMAL, 'Local DRW file exists', drwfile.url())

            # ... otherwise append DRW to container
            else:

                # Set DRW energy range
                drw.ebounds(gammalib.GEbounds(ebounds.emin(i), ebounds.emax(i)))

                # Append DRW to container
                drws.append(drw)

                # Append energy index to container
                engindex.append(i)

                # Append filename to list
                drwfiles.append(drwfile)

                # Log appending
                self._log_value(gammalib.NORMAL, 'Append DRW for computation', drwfile.url())

            # Append DRW filename in output folder
            drwnames.append(drwname)

        # If there are DRW files to compute then compute and save them now
        if len(drwfiles) > 0:

            # Write header
            self._log_header3(gammalib.NORMAL, 'Compute DRWs')

            # Compute DRWs
            if drwmethod == 'const':

                # Load DRG
                drg = gammalib.GCOMDri(drgname)

                # Get dataspace dimensions
                nphibar = drg.nphibar()
                npix    = drg.nchi() * drg.npsi()

                # Set DRWs by multiplying DRG with solid angle
                for ipix in range(npix):
                    omega = drg.map().solidangle(ipix)
                    for iphibar in range(nphibar):
                        index = ipix + iphibar*npix
                        drw   = drg[index] * omega
                        for i in range(drws.size()):
                            drws[i][index] = drw

            else:
                drws.compute_drws(obs, self._select, self['zetamin'].real(),
                                                     self['timebin'].real(),
                                                     drwmethod)

            # Phibar normalise DRWs to DREs
            for i in range(drws.size()):

                # Load DRE
                dre = gammalib.GCOMDri(drenames[engindex[i]])

                # Get dataspace dimensions
                nphibar = dre.nphibar()
                npix    = dre.nchi() * dre.npsi()

                # Phibar normalise DRW
                for iphibar in range(nphibar):
                    sum_dre = 0.0
                    sum_drw = 0.0
                    for ipix in range(npix):
                        index    = ipix + iphibar*npix
                        sum_dre += dre[index]
                        sum_drw += drws[i][index]
                    if sum_drw > 0:
                        for ipix in range(npix):
                            index           = ipix + iphibar*npix
                            drws[i][index] *= sum_dre / sum_drw

            # Save DRWs
            for i, drwfile in enumerate(drwfiles):

                # Write header
                self._log_header3(gammalib.NORMAL, 'Save DRW for %.3f - %.3f MeV' % \
                                  (ebounds.emin(i).MeV(), ebounds.emax(i).MeV()))

                # Save DRW
                drws[i].save(drwfile, True)

                # Log creation
                self._log_value(gammalib.NORMAL, 'Local DRW file created', drwfile.url())

                # Log DRW
                self._log_string(gammalib.NORMAL, str(drws[i]))

        # Return
        return drxname, drgname, drwnames, drenames

    def _generate_drb(self, drename, drgname):
        """
        Generate DRB by Phibar-normalising the DRG to the DRE

        Parameters
        ----------
        drename : str
            DRE filename
        drgname : str
            DRG filename

        Returns
        -------
        drbname : str
            DRB filename
        """ 
        # Set DRB name (always in local datastore)
        drbname = drename.replace('_dre_', '_drb_').replace(self._global_datastore, self['outfolder'].string())
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
                    index    = i + i2 * npix
                    sum_dre += dre[index]
                    sum_drb += drb[index]
                if sum_drb > 0:
                    for i in range(npix):
                        index       = i + i2 * npix
                        drb[index] *= sum_dre / sum_drb

            # Save DRB
            drb.save(drbfile, True)

            # Log creation
            self._log_value(gammalib.NORMAL, 'DRB file created', drbfile.url())

            # Log DRB
            self._log_string(gammalib.NORMAL, str(drb))

        # Return
        return drbname

    def _generate_iaq(self, ebounds):
        """
        Generate IAQ from ebounds

        Parameters
        ----------
        ebounds : `~gammalib.GEbounds'
            Energy boundaries

        Returns
        -------
        iaqname : str
            IAQ filename
        """
        # Get dataspace dimension
        dchi    = self['dchi'].real()
        dpsi    = self['dpsi'].real()
        dphibar = self['dphibar'].real()
        nphibar = self['nphibar'].integer()

        # Get calibration database and response type
        caldb    = gammalib.GCaldb('cgro', 'comptel')
        response = gammalib.toupper(self['response'].string())

        # Initialse list of IAQ names
        iaqnames = []

        # Generate one IAQ for each energy boundary
        for i in range(ebounds.size()):

            # Handle simulated response
            if response == 'SIM2' or response == 'SIM3':

                # Write header
                self._log_header3(gammalib.NORMAL,
                                  'Use CALDB IAQ for %.3f - %.3f MeV' %
                                  (ebounds.emin(i).MeV(), ebounds.emax(i).MeV()))

                # Build response name
                iaqname = '%s(%.2f-%.2f)MeV(%.0f)deg' % \
                          (response, ebounds.emin(i).MeV(), ebounds.emax(i).MeV(),
                           dphibar)

                # Raise an exception if response does not exist
                filename = caldb.filename('', '', 'IAQ', '', '', iaqname)
                if filename.is_empty():
                    msg = ('Could not find IAQ "'+iaqname+'" in calibration database. '
                           'Please select a different response type.')
                    raise RuntimeError(msg)

                # Log IAQ usage
                self._log_value(gammalib.NORMAL, 'Use CALDB IAQ', iaqname)

                # Append IAQ name
                iaqnames.append(iaqname)

            # ... otherwise handle modeled response
            else:

                # Set IAQ suffix in case that the DRI definition deviates from the
                # standard
                iaq_suffix = ''
                filler     = '_'
                if dphibar != 2.0:
                    iaq_suffix += '%s%.1fdeg' % (filler, dphibar)
                    filler      = '-'
                if nphibar != 25:
                    iaq_suffix += '%s%dbins' % (filler, nphibar)
                    filler      = '-'

                # Set IAQ filename
                iaqname = '%s/iaq_%6.6d-%6.6dkeV%s.fits' % \
                          (self['outfolder'].string(),
                           ebounds.emin(i).keV(), ebounds.emax(i).keV(),
                           iaq_suffix)
                iaqfile = gammalib.GFilename(iaqname)

                # Set global IAQ filename
                iaqname_global = '%s/iaq_%6.6d-%6.6dkeV%s.fits' % \
                                 (self._global_datastore,
                                  ebounds.emin(i).keV(), ebounds.emax(i).keV(),
                                  iaq_suffix)
                iaqfile_global = gammalib.GFilename(iaqname_global)

                # Write header
                self._log_header3(gammalib.NORMAL,
                                  'Compute IAQ for %.3f - %.3f MeV' %
                                  (ebounds.emin(i).MeV(), ebounds.emax(i).MeV()))

                # If IAQ file exists in global datastore then do nothing
                if iaqfile_global.exists():
                    iaqname = iaqname_global
                    self._log_value(gammalib.NORMAL,
                                    'Global IAQ file exists', iaqfile_global.url())

                # If IAQ file exists then do nothing
                elif iaqfile.exists():
                    self._log_value(gammalib.NORMAL,
                                    'Local IAQ file exists', iaqfile.url())

                # ... otherwise compute and save it
                else:

                    # Set Phibar maximum
                    phimax = float(int(nphibar * dphibar + 0.5))

                    # Initialise IAQ
                    iaq = gammalib.GCOMIaq(55.0, dchi, phimax, dphibar)

                    # Set IAQ energy range
                    ebds = gammalib.GEbounds(ebounds.emin(i), ebounds.emax(i))

                    # Generate continuum IAQ
                    spectrum = gammalib.GModelSpectralPlaw(1.0,-2.0,gammalib.GEnergy(1.0,'MeV'))

                    # Compute IAQ
                    iaq.set(spectrum, ebds)

                    # Save IAQ
                    iaq.save(iaqfile, True)

                    # Log creation
                    self._log_value(gammalib.NORMAL, 'Local IAQ file created', iaqfile.url())

                    # Log IAQ
                    self._log_string(gammalib.NORMAL, str(iaq))

                # Append IAQ filename
                iaqnames.append(iaqname)

        # Return
        return iaqnames


    # Public methods
    def process(self):
        """
        Process the script
        """
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

        # Set selection set
        self._set_selection_set()

        # Write header
        self._log_header1(gammalib.NORMAL, 'Bin observations')

        # Get energy boundaries
        ebounds = self._get_ebounds()

        # Create observation list XML object
        xml = gammalib.GXml()

        # Append observation list
        xml_list = xml.append('observation_list title="observation list"')

        # Loop over all input observations
        for obs in self.obs():

            # Write header
            self._log_header2(gammalib.NORMAL, self._get_obs_header(obs))

            # Log observation
            self._log_string(gammalib.NORMAL, str(obs))

            # Generate DRIs
            drxname, drgname, drwnames, drenames = self._generate_dri(obs, ebounds)

            # Skip observation if no superpackets were used
            drx = gammalib.GCOMDri(drxname)
            if drx.num_used_superpackets() == 0:
                self._log_string(gammalib.NORMAL, 'No superpackets used for '
                                 'DRX file %s. Skip observation.' % drxname)
                continue

            # Generate IAQs
            iaqnames = self._generate_iaq(ebounds)

            # Loop over all energy bins
            for i, drename in enumerate(drenames):

                # Generate DRB
                drbname = self._generate_drb(drename, drgname)

                # Get IAQ name
                iaqname = iaqnames[i]

                # Set name and ID
                name = obs.name()
                id   = '%s_%6.6d-%6.6dkeV' % (obs.id(), ebounds.emin(i).keV(),
                                              ebounds.emax(i).keV())

                # Append observation
                xml_obs = xml_list.append('observation name="%s" id="%s" instrument="COM"' % \
                                          (name, id))

                # Append DRE, DRB, DRW, DRG and DRX files
                xml_obs.append('parameter name="DRE" file="%s"' % (drename))
                xml_obs.append('parameter name="DRB" file="%s"' % (drbname))
                xml_obs.append('parameter name="DRW" file="%s"' % (drwnames[i]))
                xml_obs.append('parameter name="DRG" file="%s"' % (drgname))
                xml_obs.append('parameter name="DRX" file="%s"' % (drxname))

                # Append IAQ file
                xml_obs.append('parameter name="IAQ" value="%s"' % (iaqname))

        # Read observation container from XML object
        self.obs().clear()
        self.obs().read(xml)

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
    app = comobsbin(sys.argv)

    # Execute application
    app.execute()
