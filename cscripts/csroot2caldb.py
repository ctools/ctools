#! /usr/bin/env python
# ==========================================================================
# Generate IRFs in CALDB format from a ROOT offaxis performance file
#
# Copyright (C) 2016-2017 Juergen Knoedlseder
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
import os
import sys
import math
import copy
from datetime import datetime
import gammalib
import ctools
from cscripts import calutils

# Optional ROOT import
try:
    from ROOT import TFile
    _has_root = True
except ImportError:
    _has_root = False


# ================== #
# csroot2caldb class #
# ================== #
class csroot2caldb(ctools.cscript):
    """
    Generates IRFs in CALDB from a ROOT offaxis performance file
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'csroot2caldb'
        self._version = '1.3.0'

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query parameters
        self['infile'].filename()
        self['outdir'].string()
        self['inst'].string()
        self['id'].string()
        self['analysis'].string()
        self['zenith'].real()
        self['azimuth'].real()
        self['version'].string()
        self['psftype'].string()
        self['split'].boolean()
        self['norm1d'].boolean()
        self['rebin'].boolean()
        self['eascale'].real()
        self['bgdscale'].real()
        self['bgdoversample'].integer()
        self['bgdinfill'].boolean()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _init_metadata(self):
        """
        Initialise dictionary containing the IRF metadata

        Returns
        -------
        irf : dict
            IRF metadata dictionary
        """
        # Initialise IRF metadata dictionary
        irf = {}

        # Set calibration information
        cal_name     = 'NAME('+self['id'].string()+')'
        cal_version  = 'VERSION('+self['version'].string()+')'
        cal_cut      = 'CLASS(BEST)'
        cal_analysis = 'ANALYSIS('+self['analysis'].string()+')'
        cal_zenith   = 'ZENITH(%.3f)deg' % self['zenith'].real()
        cal_azimuth  = 'AZIMUTH(%.3f)deg' % self['azimuth'].real()
        cal_bounds   = [cal_name, cal_version, cal_cut, cal_analysis, \
                        cal_zenith, cal_azimuth]

        # Set IRF information
        irf['CAL_TEL']      = 'CTA'
        irf['CAL_INST']     = self['inst'].string().upper()
        irf['CAL_OBSID']    = self['id'].string()
        irf['CAL_DET']      = 'NONE'
        irf['CAL_FLT']      = 'NONE'
        irf['CAL_CLASS']    = 'BCF'
        irf['CAL_TYPE']     = 'DATA'
        irf['CAL_QUAL']     = 0            # 0=good, 1=bad, 2=dubious, ...
        irf['CAL_DATE']     = '14/01/30'
        irf['VAL_DATE']     = '2014-01-30'
        irf['VAL_TIME']     = '00:00:00'
        irf['REF_TIME']     = 51544.0
        irf['EA_NAME']      = 'EFF_AREA'
        irf['EA_DOC']       = 'CAL/GEN/92-019'
        irf['EA_BOUNDS']    = copy.deepcopy(cal_bounds)
        irf['EA_DESC']      = 'CTA effective area'
        irf['PSF_NAME']     = 'RPSF'
        irf['PSF_DOC']      = 'CAL/GEN/92-020'
        irf['PSF_BOUNDS']   = copy.deepcopy(cal_bounds)
        irf['PSF_DESC']     = 'CTA point spread function'
        irf['EDISP_NAME']   = 'EDISP'
        irf['EDISP_DOC']    = '???'
        irf['EDISP_BOUNDS'] = copy.deepcopy(cal_bounds)
        irf['EDISP_DESC']   = 'CTA energy dispersion'
        irf['BGD_NAME']     = 'BGD'
        irf['BGD_DOC']      = '???'
        irf['BGD_BOUNDS']   = copy.deepcopy(cal_bounds)
        irf['BGD_DESC']     = 'CTA background'

        # Return metadata
        return irf

    def _make_dirs(self, version, irf):
        """
        Generate CALDB directory structure for one observation identifier

        The structure is given by

            data/<tel>/<inst>/bcf/<obsid>

        where <tel> is "cta" and <inst> is the instrument specified in the
        CALDB constructor (the instrument may be used for different array
        configurations).

        Parameters
        ----------
        version : str
            Version string
        irf : dict
            IRF metadata dictionary

        Returns
        -------
        ds : dict
            Directory structure information
        """
        # Initialise directory structure dictionary
        ds = {}

        # Set base directory
        base_dir  = 'data'
        base_dir += '/'+irf['CAL_TEL'].lower()
        base_dir += '/'+irf['CAL_INST'].lower()

        # Set directory names
        obs_dir          = base_dir+'/bcf/'+irf['CAL_OBSID']
        ds['EA_DIR']    = obs_dir
        ds['PSF_DIR']   = obs_dir
        ds['EDISP_DIR'] = obs_dir
        ds['BGD_DIR']   = obs_dir

        # Set path to response components. If the "outdir" parameter is set
        # then prefix the value to the path of the response components
        outdir = self['outdir'].string()
        if len(outdir) > 0:
            ds['BASE_PATH']  = outdir+'/'+base_dir
            ds['EA_PATH']    = outdir+'/'+ds['EA_DIR']
            ds['PSF_PATH']   = outdir+'/'+ds['PSF_DIR']
            ds['EDISP_PATH'] = outdir+'/'+ds['EDISP_DIR']
            ds['BGD_PATH']   = outdir+'/'+ds['BGD_DIR']
        else:
            ds['BASE_PATH']  = base_dir
            ds['EA_PATH']    = ds['EA_DIR']
            ds['PSF_PATH']   = ds['PSF_DIR']
            ds['EDISP_PATH'] = ds['EDISP_DIR']
            ds['BGD_PATH']   = ds['BGD_DIR']

        # Set IRF component file names. If the "split" parameter is "yes" then
        # split the IRF components over several files.
        if self['split'].boolean():
            ds['EA_FILE']    = 'ea_'+version+'.fits'
            ds['PSF_FILE']   = 'psf_'+version+'.fits'
            ds['EDISP_FILE'] = 'edisp_'+version+'.fits'
            ds['BGD_FILE']   = 'bgd_'+version+'.fits'
        else:
            ds['EA_FILE']    = 'irf_'+version+'.fits'
            ds['PSF_FILE']   = 'irf_'+version+'.fits'
            ds['EDISP_FILE'] = 'irf_'+version+'.fits'
            ds['BGD_FILE']   = 'irf_'+version+'.fits'

        # Create directory structure
        if not os.path.isdir(ds['EA_PATH']):
            os.makedirs(ds['EA_PATH'])
        if not os.path.isdir(ds['PSF_PATH']):
            os.makedirs(ds['PSF_PATH'])
        if not os.path.isdir(ds['EDISP_PATH']):
            os.makedirs(ds['EDISP_PATH'])
        if not os.path.isdir(ds['BGD_PATH']):
            os.makedirs(ds['BGD_PATH'])

        # Return ds
        return ds

    def _open(self, irf, ds):
        """
        Open existing or create new calibration

        The actual version will put all calibrations in the same file, although
        each part of the response function will have its own logical name. We
        can thus easily modify the script to put each calibration information
        in a separate file

        Parameters
        ----------
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Open calibration database index
        ds['CIF'] = gammalib.GFits(ds['BASE_PATH']+'/caldb.indx', True)

        # If file has no CIF extension than create it now
        if not ds['CIF'].contains('CIF'):
            ds['CIF'].append(calutils.create_cif_table())
        ds['HDU_CIF'] = ds['CIF'].table('CIF')

        # Set IRF component filenames
        ea_filename    = ds['EA_PATH']+'/'+ds['EA_FILE']
        psf_filename   = ds['PSF_PATH']+'/'+ds['PSF_FILE']
        edisp_filename = ds['EDISP_PATH']+'/'+ds['EDISP_FILE']
        bgd_filename   = ds['BGD_PATH']+'/'+ds['BGD_FILE']

        # Open IRF component files
        if self['split'].boolean():
            ds['EA_FITS']    = gammalib.GFits(ea_filename, True)
            ds['PSF_FITS']   = gammalib.GFits(psf_filename, True)
            ds['EDISP_FITS'] = gammalib.GFits(edisp_filename, True)
            ds['BGD_FITS']   = gammalib.GFits(bgd_filename, True)
        else:
            ds['EA_FITS']    = gammalib.GFits(ea_filename, True)
            ds['PSF_FITS']   = ds['EA_FITS']
            ds['EDISP_FITS'] = ds['EA_FITS']
            ds['BGD_FITS']   = ds['EA_FITS']

        # Open HDUs
        ds['HDU_EA']    = self._open_hdu(ds['EA_FITS'], "EFFECTIVE AREA",
                                         irf['EA_NAME'], irf['EA_DOC'],
                                         irf)
        ds['HDU_PSF']   = self._open_hdu(ds['PSF_FITS'], "POINT SPREAD FUNCTION",
                                         irf['PSF_NAME'], irf['PSF_DOC'],
                                         irf)
        ds['HDU_EDISP'] = self._open_hdu(ds['EDISP_FITS'], "ENERGY DISPERSION",
                                         irf['EDISP_NAME'], irf['EDISP_DOC'],
                                         irf)
        ds['HDU_BGD']   = self._open_hdu(ds['BGD_FITS'], "BACKGROUND",
                                         irf['BGD_NAME'], irf['BGD_DOC'],
                                         irf)

        # Return
        return

    def _close(self, irf, ds):
        """
        Close calibration FITS files

        Parameters
        ----------
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Add information to CIF. We do this now as before we did not
        # necessarily have all the information at hand (in particular about
        # the boundaries)
        self._add_cif_info(irf, ds)

        # Save and close CIF
        ds['CIF'].save(True)
        ds['CIF'].close()

        # Close all IRF components. If the files have been split all components
        # need to be closed, otherwise only the effective area component needs
        # to be closed as representative for all components.
        ds['EA_FITS'].save(True)
        ds['EA_FITS'].close()
        if self['split'].boolean():
            ds['PSF_FITS'].save(True)
            ds['EDISP_FITS'].save(True)
            ds['BGD_FITS'].save(True)
            ds['PSF_FITS'].close()
            ds['EDISP_FITS'].close()
            ds['BGD_FITS'].close()
 
        # Return
        return

    def _open_hdu(self, fits, extname, name, doc, irf):
        """
        Open HDU

        Opens a FITS binary table with given "extname". If HDU does not exist
        in the FITS file it will be created and appended to the FITS file.

        Parameters
        ----------
        fits : `~gammalib.GFits`
            FITS file
        extname : str
            Extension name
        name : str
            Name string
        doc : str
            Document string
        irf : dict
            IRF metadata dictionary

        Returns
        -------
        table : `~gammalib.GFitsBinTable`
            FITS binary table
        """
        # Create table if it does not yet exist
        if not fits.contains(extname):

            # Create binary table
            table = gammalib.GFitsBinTable()

            # Set extension name
            table.extname(extname)

            # Set OGIP keywords
            self._set_ogip_keywords(table, doc, ['RESPONSE', name], irf)

            # Append table to FITS file
            fits.append(table)

        # Return FITS table
        return fits.table(extname)

    def _set_ogip_keywords(self, hdu, hdudoc, hduclas, irf):
        """
        Set standard OGIP keywords for extension

        Parameters
        ----------
        hdu : `~gammalig.GFitsHDU`
            Header Data Unit
        hdudoc : str
            Documentation reference string
        hduclas : list of str
            List of HDUCLAS fields
        irf : dict
            IRF dictonary
        """
        # Set UTC date of file creation
        utc = datetime.utcnow().isoformat()[:19]

        # Set keywords
        hdu.card('ORIGIN', 'IRAP', 'Name of organization making this file')
        hdu.card('DATE', utc, 'File creation date (YYYY-MM-DDThh:mm:ss UTC)')
        hdu.card('TELESCOP', irf['CAL_TEL'], 'Name of telescope')
        hdu.card('INSTRUME', irf['CAL_INST'], 'Name of instrument')
        hdu.card('DETNAM', irf['CAL_DET'], 'Name of detector')
        hdu.card('HDUCLASS', 'OGIP', 'HDU class')
        hdu.card('HDUDOC', hdudoc, 'HDU documentation')
        for i, item in enumerate(hduclas):
            key = 'HDUCLAS%d' % (i+1)
            hdu.card(key, item, 'HDU class')
        hdu.card('HDUVERS', '1.0.0', 'HDU version')

        # Return
        return

    def _add_cif_info(self, irf, ds):
        """
        Add information to CIF extension

        Parameters
        ----------
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Set list of IRF component names
        names = ['EA', 'PSF', 'EDISP', 'BGD']

        # Initialise CIF row index
        row = ds['HDU_CIF'].nrows()

        # Append rows for all components to CIF extension
        ds['HDU_CIF'].append_rows(len(names))

        # Add information for all components
        for name in names:

            # Set dictionary keys for this component
            key_dir    = '%s_DIR' % name
            key_file   = '%s_FILE' % name
            key_name   = '%s_NAME' % name
            key_desc   = '%s_DESC' % name
            key_bounds = '%s_BOUNDS' % name

            # Set generic information
            ds['HDU_CIF']['TELESCOP'][row] = irf['CAL_TEL']
            ds['HDU_CIF']['INSTRUME'][row] = irf['CAL_INST']
            ds['HDU_CIF']['DETNAM'][row]   = irf['CAL_DET']
            ds['HDU_CIF']['FILTER'][row]   = irf['CAL_FLT']
            ds['HDU_CIF']['CAL_DEV'][row]  = 'ONLINE'
            ds['HDU_CIF']['CAL_CLAS'][row] = irf['CAL_CLASS']
            ds['HDU_CIF']['CAL_DTYP'][row] = irf['CAL_TYPE']
            ds['HDU_CIF']['CAL_VSD'][row]  = irf['VAL_DATE']
            ds['HDU_CIF']['CAL_VST'][row]  = irf['VAL_TIME']
            ds['HDU_CIF']['REF_TIME'][row] = irf['REF_TIME']
            ds['HDU_CIF']['CAL_QUAL'][row] = irf['CAL_QUAL']
            ds['HDU_CIF']['CAL_DATE'][row] = irf['CAL_DATE']

            # Set component specific information
            ds['HDU_CIF']['CAL_DIR'][row]   = ds[key_dir]
            ds['HDU_CIF']['CAL_FILE'][row]  = ds[key_file]
            ds['HDU_CIF']['CAL_CNAM'][row]  = irf[key_name]
            ds['HDU_CIF']['CAL_DESC'][row]  = irf[key_desc]
            ds['HDU_CIF']['CAL_XNO'][row]   = 1
            for i in range(9):
                if i >= len(irf[key_bounds]):
                    ds['HDU_CIF']['CAL_CBD'][row,i] = 'NONE'
                else:
                    ds['HDU_CIF']['CAL_CBD'][row,i] = irf[key_bounds][i]

            # Increment row index
            row += 1

        # Return
        return

    def _set_cif_keywords(self, hdu, name, bounds, desc, irf):
        """
        Set standard CIF keywords for extension

        Parameters
        ----------
        hdu : `~gammalib.GFitsHDU`
            FITS HDU
        name : str
            Calibration name
        bounds : list of str
            Calibration boundaries
        desc : str
            Calibration description
        irf : dict
            IRF metadata dictionary
        """
        # Set standard CIF keywords
        hdu.card('CSYSNAME', 'XMA_POL', '')
        hdu.card('CCLS0001', irf['CAL_CLASS'], 'Calibration class')
        hdu.card('CDTP0001', irf['CAL_TYPE'], 'Calibration type')
        hdu.card('CCNM0001', name, 'Calibration name')

        # Set boundary keywords
        for i in range(9):
            keyname = 'CBD%d0001' % (i+1)
            if i >= len(bounds):
                value = 'NONE'
            else:
                value = bounds[i]
            hdu.card(keyname, value, 'Calibration boundary')

        # Set validity keywords
        hdu.card('CVSD0001', irf['VAL_DATE'],
                 'Calibration validity start date (UTC)')
        hdu.card('CVST0001', irf['VAL_TIME'],
                 'Calibration validity start time (UTC)')
        hdu.card('CDEC0001', desc, 'Calibration description')

        # Return
        return


    # ROOT specific private members
    def _root2caldb(self, irf, ds):
        """
        Translate ROOT to CALDB information

        Parameters
        ----------
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Open ROOT performance file
        tfile = TFile(self['infile'].filename().url())

        # Create effective area
        self._root2ea(tfile, irf, ds)

        # Create point spread function
        self._root2psf(tfile, irf, ds)

        # Create energy dispersion
        self._root2edisp(tfile, irf, ds)

        # Create background
        self._root2bgd(tfile, irf, ds)

        # Return
        return

    def _root2ea(self, tfile, irf, ds):
        """
        Translate ROOT to CALDB effective area extension
 
        The following ROOT histograms are used:
        - EffectiveAreaEtrue_offaxis -> EFFAREA

        Parameters
        ----------
        tfile : `~ROOT.TFile`
            ROOT file
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Generate effective area extension')

        # Get relevant ROOT histograms
        etrue = tfile.Get('EffectiveAreaEtrue_offaxis')

        # If requested then normalize the 2D histogram on the on-axis 1D
        # histogram. This assures that the 2D histogram has the same on-axis
        # effective area dependence as the 1D histogram.
        if self['norm1d'].boolean():
            etrue_1D = tfile.Get('EffectiveAreaEtrue')
            self._renorm_onaxis(etrue, etrue_1D)

        # If requested then rebin the effective area Etrue histogram. This
        # increases the number of bins by a factor 10.
        if self['rebin'].boolean():
            etrue.RebinX(10)
            neng    = etrue.GetXaxis().GetNbins()
            noffset = etrue.GetYaxis().GetNbins()
            for ioff in range(noffset):
                for ieng in range(neng):
                    value = etrue.GetBinContent(ieng+1,ioff+1) / 10.0
                    etrue.SetBinContent(ieng+1,ioff+1,value)

        # Get effective area multiplicator. This allows renormalising the
        # effective area.
        eascale = self['eascale'].real()

        # If renormalisation has been requested then do it now
        if eascale != 1.0:
            neng    = etrue.GetXaxis().GetNbins()
            noffset = etrue.GetYaxis().GetNbins()
            for ioff in range(noffset):
                for ieng in range(neng):
                    value = etrue.GetBinContent(ieng+1,ioff+1) * eascale
                    etrue.SetBinContent(ieng+1,ioff+1,value)

        # Write boundary keywords
        self._set_cif_keywords(ds['HDU_EA'], irf['EA_NAME'],
                               irf['EA_BOUNDS'], irf['EA_DESC'], irf)

        # Optionally write energy thresholds
        if self['emin'].is_valid():
            emin = self['emin'].real()
            ds['HDU_EA'].card('LO_THRES', emin, '[TeV] Low energy threshold')
        if self['emax'].is_valid():
            emin = self['emax'].real()
            ds['HDU_EA'].card('HI_THRES', emin, '[TeV] High energy threshold')

        # Create "EFFAREA" data column
        self._make_2D(etrue, ds['HDU_EA'], 'EFFAREA', 'm2')

        # Return
        return

    def _root2psf(self, tfile, irf, ds):
        """
        Translate ROOT to CALDB point spread function extension

        Parameters
        ----------
        tfile : `~ROOT.TFile`
            ROOT file
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # King profile PSF
        if self['psftype'].string() == 'King':
            self._root2psf_king(tfile, irf, ds)

        # ... otherwise use Gaussian profile PSF
        else:
            self._root2psf_gauss(tfile, irf, ds)

        # Return
        return

    def _root2psf_gauss(self, tfile, irf, ds):
        """
        Translate ROOT to CALDB point spread function extension

        The following ROOT histograms are used:
        - 1/(2*pi*SIGMA_1) -> SCALE
        - AngRes_offaxis -> SIGMA_1 (scaling: 1/0.8)
        - 0.0 -> AMPL_2
        - 0.0 -> SIGMA_2
        - 0.0 -> AMPL_3
        - 0.0 -> SIGMA_3

        Parameters
        ----------
        tfile : `~ROOT.TFile`
            ROOT file
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Write header
        self._log_header1(gammalib.TERSE,
                          'Generate Gaussian point spread function extension')

        # Get relevant ROOT histograms
        r68 = tfile.Get('AngRes_offaxis')

        # Extract number of bins in histogram
        neng    = r68.GetXaxis().GetNbins()
        noffset = r68.GetYaxis().GetNbins()

        # Converts 68% -> 1 sigma
        r68_to_sigma = 0.6624305
        for ioff in range(noffset):
            for ieng in range(neng):
                sigma = r68.GetBinContent(ieng+1,ioff+1) * r68_to_sigma
                r68.SetBinContent(ieng+1,ioff+1,sigma)

        # Compute scale histogram
        scale = r68.Clone()
        for ioff in range(noffset):
            for ieng in range(neng):
                integral = 2.0 * math.pi * r68.GetBinContent(ieng+1,ioff+1)
                if integral > 0.0:
                    value = 1.0 / integral
                else:
                    value = 0.0
                scale.SetBinContent(ieng+1,ioff+1,value)

        # Set zero histogram
        zero = r68.Clone()
        for ioff in range(noffset):
            for ieng in range(neng):
                zero.SetBinContent(ieng+1,ioff+1,0.0)

        # Set boundaries
        irf['PSF_BOUNDS'].append('PSF(GAUSS)')

        # Write boundary keywords
        self._set_cif_keywords(ds['HDU_PSF'], irf['PSF_NAME'],
                               irf['PSF_BOUNDS'], irf['PSF_DESC'], irf)

        # Create "SCALE" data column
        self._make_2D(scale, ds['HDU_PSF'], 'SCALE', '')

        # Create "SIGMA_1" data column
        self._make_2D(r68, ds['HDU_PSF'], 'SIGMA_1', 'deg')

        # Create "AMPL_2" data column
        self._make_2D(zero, ds['HDU_PSF'], 'AMPL_2', '')

        # Create "SIGMA_2" data column
        self._make_2D(zero, ds['HDU_PSF'], 'SIGMA_2', 'deg')

        # Create "AMPL_3" data column
        self._make_2D(zero, ds['HDU_PSF'], 'AMPL_3', '')

        # Create "SIGMA_3" data column
        self._make_2D(zero, ds['HDU_PSF'], 'SIGMA_3', 'deg')

        # Return
        return

    def _root2psf_king(self, tfile, irf, ds):
        """
        Translate ROOT to CALDB point spread function extension

        The following ROOT histograms are used:
        - AngRes_offaxis
        - AngRes80_offaxis

        Parameters
        ----------
        tfile : `~ROOT.TFile`
            ROOT file
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Write header
        self._log_header1(gammalib.TERSE,
                          'Generate King point spread function extension')

        # Get relevant ROOT histograms
        r68 = tfile.Get('AngRes_offaxis')
        r80 = tfile.Get('AngRes80_offaxis')

        # Extract number of bins in histogram
        neng    = r68.GetXaxis().GetNbins()
        noffset = r68.GetYaxis().GetNbins()

        # Initialise parameter maps by cloning the r68 2D map
        gamma2D = r68.Clone()
        sigma2D = r68.Clone()

        # Compute gamma and sigma values
        for ioff in range(noffset):

            # Initialise last results
            last_gamma = 0.0
            last_sigma = 0.0

            # Loop over all energies
            for ieng in range(neng):

                # Extract radii
                r_68 = r68.GetBinContent(ieng+1,ioff+1)
                r_80 = r80.GetBinContent(ieng+1,ioff+1)

                # Initialise results
                gamma = 0.0
                sigma = 0.0

                # Continue only if both radii are positive
                if r_68 > 0 and r_80 > 0:

                    # Derive constants for equation to solve
                    a = 1.0 - 0.68
                    b = 1.0 - 0.80
                    c = r_68*r_68/(r_80*r_80)

                    # Solve equation (a^x-1)/(b^x-1)=c for x using secant
                    # method. Stop when we are better than 1e-6.
                    x1   = -0.5
                    x2   = -1.0
                    f1   = (math.pow(a,x1) - 1.0)/(math.pow(b,x1) - 1.0) - c
                    f2   = (math.pow(a,x2) - 1.0)/(math.pow(b,x2) - 1.0) - c
                    while True:
                        x     = x1 - f1 * (x1-x2)/(f1-f2)
                        f     = (math.pow(a,x) - 1.0)/(math.pow(b,x) - 1.0) - c
                        if abs(f) < 1.0e-6:
                            break
                        else:
                            f2 = f1
                            x2 = x1
                            f1 = f
                            x1 = x

                    # Compute gamma.
                    if x < 0.0:
                        gamma = 1.0 - 1.0/x
                    else:
                        gamma = 1.0

                    # Compute sigma
                    denominator = 2.0 * gamma * (math.pow(a, x) - 1.0)
                    if denominator > 0.0:
                        sigma = r_68 * math.sqrt(1.0/denominator)
                    else:
                        denominator = 2.0 * gamma * (math.pow(b, x) - 1.0)
                        if denominator > 0.0:
                            sigma = r_80 * math.sqrt(1.0/denominator)
                        else:
                            gamma = 0.0
                            sigma = 0.0

                    # Handle special case that no results were found.
                    # This takes care of pixels that are ill defined
                    # in the MC file.
                    if gamma == 0.0 and sigma == 0.0:
                        gamma  = last_gamma
                        sigma  = last_sigma
                        status = ' (use last)'
                    else:
                        status = ''

                    # Log results
                    self._log_value(gammalib.EXPLICIT,
                                    'ieng=%d ioff=%d%s' % (ieng,ioff,status),
                                    'r68=%f r80=%f gamma=%f sigma=%f' %
                                    (r_68, r_80, gamma, sigma))
    
                    # Store surrent result as last result
                    last_gamma = gamma
                    last_sigma = sigma

                # Set bin contents
                gamma2D.SetBinContent(ieng+1,ioff+1,gamma)
                sigma2D.SetBinContent(ieng+1,ioff+1,sigma)

        # Set boundaries
        irf['PSF_BOUNDS'].append('PSF(KING)')

        # Write boundary keywords
        self._set_cif_keywords(ds['HDU_PSF'], irf['PSF_NAME'],
                               irf['PSF_BOUNDS'], irf['PSF_DESC'], irf)

        # Create "GAMMA" data column
        self._make_2D(gamma2D, ds['HDU_PSF'], 'GAMMA', '')

        # Create "SIGMA" data column
        self._make_2D(sigma2D, ds['HDU_PSF'], 'SIGMA', 'deg')

        # Return
        return

    def _root2edisp(self, tfile, irf, ds):
        """
        Translate ROOT to CALDB energy dispersion extension

        The following ROOT histograms are used:
        - EestOverEtrue_offaxis  -> MATRIX

        Parameters
        ----------
        tfile : `~ROOT.TFile`
            ROOT file
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Generate energy dispersion extension')

        # Get relevant ROOT histograms
        matrix = tfile.Get('EestOverEtrue_offaxis')

        # Write boundary keywords
        self._set_cif_keywords(ds['HDU_EDISP'], irf['EDISP_NAME'],
                               irf['EDISP_BOUNDS'], irf['EDISP_DESC'], irf)

        # Create "MATRIX" data column
        self._make_3D_migra(matrix, ds['HDU_EDISP'], 'MATRIX', '')

        # Return
        return

    def _root2bgd(self, tfile, irf, ds):
        """
        Translate ROOT to CALDB background extension.

        The following ROOT histograms are used:
        - BGRatePerSqDeg_offaxis -> BGD

        Parameters
        ----------
        tfile : `~ROOT.TFile`
            ROOT file
        irf : dict
            IRF metadata dictionary
        ds : dict
            Directory structure dictionary
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Generate 3D background extension')

        # Get relevant ROOT histograms
        array = tfile.Get('BGRatePerSqDeg_offaxis')

        # Replace 2D histogram values by power law extrapolation
        self._plaw_replace(array, self['bgdethres'].real())

        # If requested then fill-in empty values in 2D histogram
        if self['bgdinfill'].boolean():
            self._infill_bgd(array)

        # If requested then normalize the 2D histogram on the on-axis 1D
        # histogram. This assures that the 2D histogram has the same on-axis
        # effective area dependence as the 1D histogram.
        if self['norm1d'].boolean():
            array_1D = tfile.Get('BGRatePerSqDeg')
            self._renorm_onaxis(array, array_1D)

        # Write boundary keywords
        self._set_cif_keywords(ds['HDU_BGD'], irf['BGD_NAME'],
                               irf['BGD_BOUNDS'], irf['BGD_DESC'], irf)

        # Create "BGD" data column
        self._make_3D(array, ds['HDU_BGD'], 'BGD', '1/s/MeV/sr')

        # Return
        return

    def _append_column_axis(self, hdu, name, unit, axis, log=False):
        """
        Append column of ROOT axis values to HDU

        Parameters
        ----------
        hdu : `~gammalib.GFitsHDU`
            FITS HDU
        name : str
            Column name
        unit : str
            Column unit
        axis : `~ROOT.TAxis`
            ROOT histogram axis
        log : bool, optional
            Axis is logarithmic
        """
        # Continue only if columns does not yet exist
        if not hdu.contains(name):

            # Write header
            self._log_header3(gammalib.TERSE, 'Append axis column "%s"' % name)

            # Get number of axis bins
            nbins = axis.GetNbins()

            # Append column and set column unit
            hdu.append(gammalib.GFitsTableFloatCol(name, 1, nbins))
            hdu[name].unit(unit)

            # Do we have a LO axis? If not we have a HI axis
            if name[-2:] == 'LO':
                low = True
            else:
                low = False

            # Fill column
            for i in range(nbins):
                if low:
                    value = axis.GetBinLowEdge(i+1)
                else:
                    value = axis.GetBinUpEdge(i+1)
                if log:
                    value = pow(10.0, value)
                hdu[name][0,i] = value

            # Log values
            self._log_value(gammalib.NORMAL, 'Number of axis bins', nbins)
            self._log_value(gammalib.NORMAL, 'Unit', unit)

        # Return
        return

    def _append_column_values(self, hdu, name, unit, values):
        """
        Append column of values to HDU

        Parameters
        ----------
        hdu : `~gammalib.GFitsHDU`
            FITS HDU
        name : str
            Column name
        unit : str
            Column unit
        values : list of float
            Axis values
        """
        # Continue only if columns does not yet exist
        if not hdu.contains(name):

            # Write header
            self._log_header3(gammalib.TERSE, 'Append value column "%s"' % name)

            # Get number of values
            nbins = len(values)

            # Append column and set column unit
            hdu.append(gammalib.GFitsTableFloatCol(name, 1, nbins))
            hdu[name].unit(unit)

            # Fill column
            for i in range(nbins):
                hdu[name][0,i] = values[i]

            # Log values
            self._log_value(gammalib.NORMAL, 'Number of values', nbins)
            self._log_value(gammalib.NORMAL, 'Unit', unit)

        # Return
        return

    def _renorm_onaxis(self, hist2D, hist1D):
        """
        Renormalise 2D histogram (energy,offset) on 1D histogram (energy)

        This method makes sure that a 2D histogram has the same onaxis values
        as the corresponding 1D histogram.

        Parameters
        ----------
        hist2D : `~ROOT.TH2F`
            ROOT 2D histogram
        hist1D : `~ROOT.TH1F`
            ROOT 1D histogram
        """
        # Get 2D dimensions
        neng    = hist2D.GetXaxis().GetNbins()
        noffset = hist2D.GetYaxis().GetNbins()

        # Continue only if the 1D and 2D histograms have the same number
        # of energy bins
        if neng != hist1D.GetXaxis().GetNbins():

            # Log if energy binning differs
            self._log_value(gammalib.TERSE, 'On-axis normalisation',
                            'Impossible since energy binning of 1D histogram '
                            '"%s" does not match that of 2D histogram "%s".' %
                            (hist1D.GetName(), hist2D.GetName()))
        else:
            # Log on-axis normalisation
            self._log_value(gammalib.TERSE, 'On-axis normalisation',
                            'On-axis values of 1D histogram "%s" imposed on '
                            '2D histogram "%s".' %
                            (hist1D.GetName(), hist2D.GetName()))

            # Get a copy of 2D histogram
            hist2D_copy = hist2D.Clone()

            # Divide 2D histogram by onaxis values to remove the energy
            # dependence
            for ieng in range(neng):
                onaxis = hist1D.GetBinContent(ieng+1)
                if onaxis > 0.0:
                    for ioff in range(noffset):
                        value = hist2D.GetBinContent(ieng+1,ioff+1) / onaxis
                        hist2D.SetBinContent(ieng+1,ioff+1,value)

            # Smooth in energy direction to reduce statistical fluctuations
            for ieng in range(neng):
                for ioff in range(noffset):
                    if ieng == 0:
                        value = (2.0 * hist2D.GetBinContent(ieng+1,ioff+1) +
                                 1.0 * hist2D.GetBinContent(ieng+2,ioff+1)) / 3.0
                    elif ieng == neng-1:
                        value = (2.0 * hist2D.GetBinContent(ieng+1,ioff+1) +
                                 1.0 * hist2D.GetBinContent(ieng+0,ioff+1)) / 3.0
                    else:
                        value = (2.0 * hist2D.GetBinContent(ieng+1,ioff+1) +
                                 1.0 * hist2D.GetBinContent(ieng+0,ioff+1) +
                                 1.0 * hist2D.GetBinContent(ieng+2,ioff+1)) / 4.0
                    hist2D_copy.SetBinContent(ieng+1,ioff+1,value)

            # Normalize now for an onaxis value of 1
            for ieng in range(neng):
                onaxis = hist2D_copy.GetBinContent(ieng+1,1)
                if onaxis > 0.0:
                    for ioff in range(noffset):
                        value = hist2D_copy.GetBinContent(ieng+1,ioff+1) / onaxis
                        hist2D_copy.SetBinContent(ieng+1,ioff+1,value)

            # Put back the energy dependence
            for ieng in range(neng):
                onaxis = hist1D.GetBinContent(ieng+1)
                for ioff in range(noffset):
                    value = hist2D_copy.GetBinContent(ieng+1,ioff+1) * onaxis
                    hist2D.SetBinContent(ieng+1,ioff+1,value)

        # Return
        return

    def _plaw_replace(self, hist2D, ethres):
        """
        Replace background rate values by power law values

        Parameters
        ----------
        hist2D : `~ROOT.TH2F`
            ROOT 2D histogram
        ethres : float
            Energy threshold
        """
        # Extract energy and offset angle vectors
        energies = hist2D.GetXaxis()
        neng     = energies.GetNbins()
        offsets  = hist2D.GetYaxis()
        noffset  = offsets.GetNbins()

        # Continue only if threshold is low enough
        if ethres < math.pow(10.0, energies.GetBinUpEdge(neng)):

            # Determine mean energy values
            engs = []
            for ieng in range(neng):
                energy = 0.5 * (math.pow(10.0, energies.GetBinLowEdge(ieng+1)) +
                                math.pow(10.0, energies.GetBinUpEdge(ieng+1)))
                engs.append(energy)

            # Loop over all offaxis angles
            for ioff in range(noffset):

                # Determine offset angle (for logging)
                offset = 0.5 * (offsets.GetBinLowEdge(ioff+1) +
                                offsets.GetBinUpEdge(ioff+1))

                # Initialise vectors
                axis  = []
                rates = []

                # Extract values
                for ieng in range(neng):
                    energy = engs[ieng]
                    rate   = hist2D.GetBinContent(ieng+1,ioff+1)
                    axis.append(energy)
                    rates.append(rate)

                # Determine minimum energy (ignoring the maximum energy as
                # we do not need it)
                emin = self._plaw_energy_range(axis, rates)[0]

                # Get power law coefficients
                coeff = self._plaw_coefficients(axis, rates, emin, ethres)

                # Replace histogram values by power law values if coefficients
                # are valid
                if coeff['m'] != 0.0 and coeff['t'] != 0.0:
                    for ieng in range(neng):
                        energy = engs[ieng]
                        if energy > ethres:
                            value = self._plaw_value(coeff, energy)
                            hist2D.SetBinContent(ieng+1,ioff+1,value)

                            # Log power law replacement
                            self._log_value(gammalib.NORMAL,
                                            'Power law replacement',
                                            '%e cts/s/deg2 (E=%8.4f TeV (%d) '
                                            'Off=%.2f deg (%d)' %
                                            (value,engs[ieng],ieng,offset,ioff))

        # Return
        return

    def _plaw_energy_range(self, energies, rates):
        """
        Determine the energy range over which the power law is valid

        Parameters
        ----------
        energies : list of float
            Energy values (TeV)
        rates : list of float
            Background rate values
        """
        # Find energy at which background rate takes a maximum
        min_energy = 0.0
        max_rate   = 0.0
        for ieng, energy in enumerate(energies):
            if rates[ieng] > max_rate:
                max_rate   = rates[ieng]
                min_energy = energy

        # Find last data point
        max_energy = 0.0
        for ieng in range(len(energies)-1,-1,-1):
            if rates[ieng] > 0.0:
                max_energy = energies[ieng]
                break

        # Return
        return (min_energy, max_energy)

    def _plaw_coefficients(self, energies, rates, emin, emax):
        """
        Determines the power law coefficients

        Parameters
        ----------
        energies : list of float
            Energy values (TeV)
        rates : list of float
            Background rate values
        emin : float
            Minimum energy (TeV)
        emax : float
            Maximum energy (TeV)
        """
        # Determine nearest energy indices
        i_emin     = 0
        i_emax     = 0
        delta_emin = 1.0e30
        delta_emax = 1.0e30
        for ieng, energy in enumerate(energies):
            if abs(energy-emin) < delta_emin:
                delta_emin = abs(energy-emin)
                i_emin     = ieng
            if abs(energy-emax) < delta_emax:
                delta_emax = abs(energy-emax)
                i_emax     = ieng

        # Determine coefficients of power law
        while rates[i_emax] == 0.0 and i_emax > i_emin:
            i_emax -= 1
        if rates[i_emax] > 0.0:
            x1 = math.log10(energies[i_emin])
            x2 = math.log10(energies[i_emax])
            y1 = math.log10(rates[i_emin])
            y2 = math.log10(rates[i_emax])
            m  = (y2-y1)/(x2-x1)
            t  = y1 - m * x1
        else:
            m = 0.0
            t = 0.0

        # Return coefficients
        return {'m': m, 't': t}

    def _plaw_value(self, coeff, energy):
        """
        Determine the power law value

        Parameters
        ----------
        coeff : dict
            Dictionary of coefficients
        energy : float
            Energy value (TeV)
        """
        # Compute value
        value = math.pow(10.0, coeff['m'] * math.log10(energy) + coeff['t'])

        # Return value
        return value

    def _infill_bgd(self, hist2D):
        """
        Fill all empty bins in 2D histogram by preceeding values in energy

        Parameters
        ----------
        hist2D : `~ROOT.TH2F`
            ROOT 2D histogram
        """
        # Extract energy and offset angle vectors
        energies = hist2D.GetXaxis()
        neng     = energies.GetNbins()
        offsets  = hist2D.GetYaxis()
        noffset  = offsets.GetNbins()

        # Determine mean energy values (for logging)
        engs = []
        for ieng in range(neng):
            energy = 0.5 * (math.pow(10.0, energies.GetBinLowEdge(ieng+1)) +
                            math.pow(10.0, energies.GetBinUpEdge(ieng+1)))
            engs.append(energy)

        # Loop over all offaxis angles
        for ioff in range(noffset):

            # Determine offset angle (for logging)
            offset = 0.5 * (offsets.GetBinLowEdge(ioff+1) +
                            offsets.GetBinUpEdge(ioff+1))

            # Initialise indices and value
            i_start = -1
            i_stop  = -1
            value   = 0.0

            # First start and stop indices for infill
            for ieng in range(neng):
                if hist2D.GetBinContent(ieng+1,ioff+1) > 0.0:
                    i_start = ieng
                    value   = hist2D.GetBinContent(ieng+1,ioff+1)
                    break
            for ieng in range(neng-1, -1, -1):
                if hist2D.GetBinContent(ieng+1,ioff+1) > 0.0:
                    i_stop = ieng
                    break

            # If indices are valid then fill in background rates
            if i_start > -1 and i_stop > -1:
                for ieng in range(i_start, i_stop+1):
                    if hist2D.GetBinContent(ieng+1,ioff+1) == 0.0:
                        hist2D.SetBinContent(ieng+1,ioff+1,value)

                        # Log background infill
                        self._log_value(gammalib.NORMAL,
                                        'Background infill',
                                        '%e cts/s/deg2 (E=%8.4f TeV (%d) '
                                        'Off=%.2f deg (%d)' %
                                        (value,engs[ieng],ieng,offset,ioff))

                    else:
                        value = hist2D.GetBinContent(ieng+1,ioff+1)

        # Return
        return

    def _make_2D(self, array, hdu, name, unit, scale=1.0):
        """
        Make 2D data column as function of energy and offset angle

        If the HDU has already the energy and offset angle columns, this method
        will simply add another data column. If name==None, the method will not
        append any data column.

        Parameters
        ----------
        array : `~ROOT.TH2F`
            ROOT 2D histogram
        hdu : `~gammalib.GFitsHDU`
            FITS HDU
        name : str
            Data column name
        unit : str
            Data unit
        scale : float, optional
            Scaling factor for histogram values
        """
        # Write header
        self._log_header3(gammalib.TERSE, 'Make 2D data column "%s"' % name)

        # Extract energy and offset angle vectors
        energies = array.GetXaxis()
        offsets  = array.GetYaxis()
        neng     = energies.GetNbins()
        noffset  = offsets.GetNbins()

        # Log parameters
        self._log_value(gammalib.NORMAL, 'Number of energies', neng)
        self._log_value(gammalib.NORMAL, 'Number of offsets', noffset)

        # Append axis columns to HDU
        self._append_column_axis(hdu, 'ENERG_LO', 'TeV', energies, log=True)
        self._append_column_axis(hdu, 'ENERG_HI', 'TeV', energies, log=True)
        self._append_column_axis(hdu, 'THETA_LO', 'deg', offsets)
        self._append_column_axis(hdu, 'THETA_HI', 'deg', offsets)

        # Append array column
        if name != None and not hdu.contains(name):
            hdu.append(gammalib.GFitsTableFloatCol(name, 1, neng*noffset))
            hdu[name].unit(unit)
            hdu[name].dim([neng, noffset])
            for ioff in range(noffset):
                for ieng in range(neng):
                    index = ieng + ioff * neng
                    value = array.GetBinContent(ieng+1,ioff+1)
                    hdu[name][0,index] = value * scale

        # Collect boundary information
        bd_eng = 'ENERG(%.4f-%.2f)TeV' % \
                 (pow(10.0, energies.GetBinLowEdge(1)), \
                  pow(10.0, energies.GetBinUpEdge(neng)))
        bd_off = 'THETA(%.2f-%.2f)deg' % \
                 (offsets.GetBinLowEdge(1), \
                  offsets.GetBinUpEdge(noffset))
        bd_phi = 'PHI(0-360)deg'
        bounds = [bd_eng, bd_off, bd_phi]

        # Return boundary information
        return bounds

    def _make_3D(self, array, hdu, name, unit):
        """
        Make 3D data column as function of DETX, DETY and energy

        If the HDU has already the energy and offset angle columns, this method
        will simply add another data column. If name==None, the method will not
        append any data column.

        Parameters
        ----------
        array : `~ROOT.TH2F`
            ROOT 2D histogram
        hdu : `~gammalib.GFitsHDU`
            FITS HDU
        name : str
            Data column name
        unit : str
            Data unit
        """
        # Write header
        self._log_header3(gammalib.TERSE, 'Make 3D data column "%s"' % name)

        # Get User parameters
        scale      = self['bgdscale'].real()
        oversample = self['bgdoversample'].integer()

        # Set constants
        deg2sr = 0.01745329*0.01745329

        # Extract energy and offset angle vectors
        energies  = array.GetXaxis()
        neng      = energies.GetNbins()
        offsets   = array.GetYaxis()
        noffset   = offsets.GetNbins()
        theta_max = offsets.GetBinUpEdge(noffset) # Maximum theta
        ewidth    = [] # in MeV
        for ieng in range(neng):
            ewidth.append(pow(10.0, energies.GetBinUpEdge(ieng+1) +6.0) -
                          pow(10.0, energies.GetBinLowEdge(ieng+1)+6.0))

        # Log parameters
        self._log_value(gammalib.NORMAL, 'Scale', scale)
        self._log_value(gammalib.NORMAL, 'Oversample', oversample)
        self._log_value(gammalib.NORMAL, 'Number of energies', neng)
        self._log_value(gammalib.NORMAL, 'Number of offsets', noffset)
        self._log_value(gammalib.NORMAL, 'Maximum offset', theta_max)

        # Build DETX and DETY axes
        ndet     = array.GetYaxis().GetNbins()
        ndets    = 2*ndet*oversample
        det_max  = array.GetYaxis().GetBinUpEdge(ndet)
        det_bin  = det_max/float(ndet*oversample)
        dets_lo  = []
        dets_hi  = []
        dets2    = []
        for i in range(ndets):
            det_lo  = -det_max + i*det_bin
            det_hi  = det_lo + det_bin
            det_val = 0.5*(det_lo + det_hi)
            dets_lo.append(det_lo)
            dets_hi.append(det_hi)
            dets2.append(det_val*det_val)

        # Log parameters
        self._log_value(gammalib.NORMAL, 'Number of DETXY bins', ndets)
        self._log_value(gammalib.NORMAL, 'Maximum DETXY', det_max)
        self._log_value(gammalib.NORMAL, 'DETXY binsize', det_bin)

        # Append DETX_LO, DETX_HI, DETY_LO and DETY_HI columns
        self._append_column_values(hdu, 'DETX_LO', 'deg', dets_lo)
        self._append_column_values(hdu, 'DETX_HI', 'deg', dets_hi)
        self._append_column_values(hdu, 'DETY_LO', 'deg', dets_lo)
        self._append_column_values(hdu, 'DETY_HI', 'deg', dets_hi)

        # Append ENERG_LO and ENERG_HI columns
        self._append_column_axis(hdu, 'ENERG_LO', 'TeV', energies, log=True)
        self._append_column_axis(hdu, 'ENERG_HI', 'TeV', energies, log=True)

        # Append array column
        if name != None and not hdu.contains(name):
            hdu.append(gammalib.GFitsTableFloatCol(name, 1, ndets*ndets*neng))
            hdu[name].unit(unit)
            hdu[name].dim([ndets,ndets,neng])
            for ix in range(ndets):
                for iy in range(ndets):
                    for ieng in range(neng):
                        index = ix + (iy + ieng * ndets) * ndets
                        theta = math.sqrt(dets2[ix] + dets2[iy])
                        if theta < theta_max:
                            binsq = array.Interpolate(energies.GetBinCenter(ieng+1),
                                                      theta)
                            value = binsq / deg2sr / (ewidth[ieng])
                            hdu[name][0,index] = value * scale

        # Collect boundary information
        bd_detx = 'DETX(%.2f-%.2f)deg' % (-det_max, det_max)
        bd_dety = 'DETY(%.2f-%.2f)deg' % (-det_max, det_max)
        bd_eng  = 'ENERG(%.4f-%.2f)TeV' % \
                  (pow(10.0, energies.GetBinLowEdge(1)),
                   pow(10.0, energies.GetBinUpEdge(neng)))
        bounds  = [bd_detx, bd_dety, bd_eng]

        # Return boundary information
        return bounds

    def _make_3D_migra(self, array, hdu, name, unit, scale=1.0):
        """
        Make 3D data column as function of ETRUE, MIGRA and THETA

        If the HDU has already the energy and offset angle columns, this method
        will simply add another data column. If name==None, the method will not
        append any data column.

        Parameters
        ----------
        array : `~ROOT.TH3F`
            ROOT 2D histogram
        hdu : `~gammalib.GFitsHDU`
            FITS HDU
        name : str
            Data column name
        unit : str
            Data unit
        scale : float, optional
            Scaling factor for histogram values
        """
        # Write header
        self._log_header3(gammalib.TERSE,
                          'Make 3D migration matrix data column "%s"' % name)

        # Extract Etrue, Eobs/Etrue and offset angle vectors
        etrue   = array.GetXaxis()
        netrue  = etrue.GetNbins()
        migra   = array.GetYaxis()
        nmigra  = migra.GetNbins()
        offsets = array.GetZaxis()
        noffset = offsets.GetNbins()
        ewidth  = [] # in MeV
        for ieng in range(netrue):
            ewidth.append(pow(10.0, etrue.GetBinUpEdge(ieng+1) +6.0) - \
                          pow(10.0, etrue.GetBinLowEdge(ieng+1)+6.0))

        # Log parameters
        self._log_value(gammalib.NORMAL, 'Number of energies', netrue)
        self._log_value(gammalib.NORMAL, 'Number of migrations', nmigra)
        self._log_value(gammalib.NORMAL, 'Number of offsets', noffset)

        # Append axis columns to HDU
        self._append_column_axis(hdu, 'ETRUE_LO', 'TeV', etrue, log=True)
        self._append_column_axis(hdu, 'ETRUE_HI', 'TeV', etrue, log=True)
        self._append_column_axis(hdu, 'MIGRA_LO', '', migra)
        self._append_column_axis(hdu, 'MIGRA_HI', '', migra)
        self._append_column_axis(hdu, 'THETA_LO', 'deg', offsets)
        self._append_column_axis(hdu, 'THETA_HI', 'deg', offsets)

        # Append migration matrix to HDU
        if name != None and not hdu.contains(name):
            hdu.append(gammalib.GFitsTableFloatCol(name, 1, netrue*nmigra*noffset))
            hdu[name].unit(unit)
            hdu[name].dim([netrue, nmigra, noffset])
            for ioff in range(noffset):
                for imigra in range(nmigra):
                    for ieng in range(netrue):
                        index = ieng + (imigra + ioff * nmigra) * netrue
                        value = array.GetBinContent(ieng+1,imigra+1,ioff+1)
                        hdu[name][0,index] = value * scale

        # Collect boundary information
        bd_eng   = 'ETRUE(%.4f-%.2f)TeV' % \
                   (pow(10.0, etrue.GetBinLowEdge(1)),
                    pow(10.0, etrue.GetBinUpEdge(netrue)))
        bd_migra = 'MIGRA(%.3f-%.3f)' % \
                   (migra.GetBinLowEdge(1),
                    migra.GetBinUpEdge(nmigra))
        bd_off   = 'THETA(%.2f-%.2f)deg' % \
                   (offsets.GetBinLowEdge(1),
                    offsets.GetBinUpEdge(noffset))
        bd_phi   = 'PHI(0-360)deg'
        bounds   = [bd_eng, bd_migra, bd_off, bd_phi]

        # Return boundary information
        return bounds


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Warn if ROOT module is missing
        if not _has_root:
            gammalib.warning('csroot2caldb', 'ROOT Python module not present, '
                             'script will create empty IRF.')

        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Initialise IRF metadata
        irf = self._init_metadata()

        # Create directory structure
        ds = self._make_dirs('file', irf)

        # Open calibration files
        self._open(irf, ds)

        # Translate ROOT to CALDB information
        if _has_root:
            self._root2caldb(irf, ds)

        # Close calibration files
        self._close(irf, ds)

        # Return
        return

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Run the script
        self.run()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csroot2caldb(sys.argv)

    # Execute application
    app.execute()
