#! /usr/bin/env python
# ==========================================================================
# Generate IRFs in CALDB format from a ROOT offaxis performance file
#
# Copyright (C) 2016 Juergen Knoedlseder
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
import gammalib
import ctools
from datetime import datetime

# Optional ROOT import
try:
    from ROOT import TFile, TH2F, TH3F
    _has_root = True
except:
    _has_root = False


# ================== #
# csroot2caldb class #
# ================== #
class csroot2caldb(ctools.cscript):
    """
    Generates IRFs in CALDB from a ROOT offaxis performance file.
    """

    # Constructor
    def __init__(self, *argv):
        """
        CALDB constructor.
        
        """
        # Set name
        self._name    = "csroot2caldb"
        self._version = "1.0.0"

        # Initialise some members
        self._cif     = None
        self._irf     = None
        self._ea      = None
        self._psf     = None
        self._edisp   = None
        self._bgd     = None
        self._hdu_cif = None

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile.
        """
        # Query parameters
        self["infile"].filename()
        self["outdir"].string()
        self["inst"].string()
        self["id"].string()
        self["analysis"].string()
        self["version"].string()
        self["psftype"].string()
        self["rebin"].boolean()
        self["eascale"].real()
        self["bgdscale"].real()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Return
        return

    def _initialise(self):

        # Set CALDB information
        self._cal_tel      = "CTA"
        self._cal_inst     = self["inst"].string().upper()
        self._cal_obsid    = self["id"].string()
        self._cal_det      = "NONE"
        self._cal_flt      = "NONE"
        self._cal_class    = "BCF"
        self._cal_type     = "DATA"
        self._cal_qual     = 0            # 0=good, 1=bad, 2=dubious, ...
        self._cal_date     = "14/01/30"
        self._val_date     = "2014-01-30"
        self._val_time     = "00:00:00"
        self._ref_time     = 51544.0
        self._cal_name     = "NAME("+self["id"].string()+")"
        self._cal_version  = "VERSION("+self["version"].string()+")"
        self._cal_cut      = "CLASS(BEST)"
        self._cal_analysis = "ANALYSIS("+self["analysis"].string()+")"
        #
        self._ea_name      = "EFF_AREA"
        self._ea_doc       = "CAL/GEN/92-019"
        self._ea_bounds    = [self._cal_name, self._cal_version,
                              self._cal_cut, self._cal_analysis]
        self._ea_desc      = "CTA effective area"
        #
        self._psf_name     = "RPSF"
        self._psf_doc      = "CAL/GEN/92-020"
        self._psf_bounds   = [self._cal_name, self._cal_version,
                              self._cal_cut, self._cal_analysis]
        self._psf_desc     = "CTA point spread function"
        #
        self._edisp_name   = "EDISP"
        self._edisp_doc    = "???"
        self._edisp_bounds = [self._cal_name, self._cal_version,
                              self._cal_cut, self._cal_analysis]
        self._edisp_desc   = "CTA energy dispersion"
        #
        self._bgd_name     = "BGD"
        self._bgd_doc      = "???"
        self._bgd_bounds   = [self._cal_name, self._cal_version,
                              self._cal_cut, self._cal_analysis]
        self._bgd_desc     = "CTA background"

        # Create directory structure
        self._make_dirs()

        # Return
        return

    def _make_dirs(self):
        """
        Generate CALDB directory structure for one observation identifier.
        The structure is given by

            data/<tel>/<inst>/bcf/<obsid>

        where <tel> is "cta" and <inst> is the instrument specified in the
        CALDB constructor (the instrument may be used for different array
        configurations).
        """
        # Set base directory
        base_dir  = "data"
        base_dir += "/"+self._cal_tel.lower()
        base_dir += "/"+self._cal_inst.lower()

        # Set directory names
        obs_dir         = base_dir+"/bcf/"+self._cal_obsid
        self._ea_dir    = obs_dir
        self._psf_dir   = obs_dir
        self._edisp_dir = obs_dir
        self._bgd_dir   = obs_dir

        # Optionally prefix path
        outdir = self["outdir"].string()
        if len(outdir) > 0:
            self._base_path  = outdir+"/"+base_dir
            self._ea_path    = outdir+"/"+self._ea_dir
            self._psf_path   = outdir+"/"+self._psf_dir
            self._edisp_path = outdir+"/"+self._edisp_dir
            self._bgd_path   = outdir+"/"+self._bgd_dir
        else:
            self._base_path  = base_dir
            self._ea_path    = self._ea_dir
            self._psf_path   = self._psf_dir
            self._edisp_path = self._edisp_dir
            self._bgd_path   = self._bgd_dir

        # Create directory structure
        if not os.path.isdir(self._ea_path):
            os.makedirs(self._ea_path)
        if not os.path.isdir(self._psf_path):
            os.makedirs(self._psf_path)
        if not os.path.isdir(self._edisp_path):
            os.makedirs(self._edisp_path)
        if not os.path.isdir(self._bgd_path):
            os.makedirs(self._bgd_path)

        # Return
        return

    def _open(self, version, split=False, clobber=True):
        """
        Open existing or create new calibration. The actual version will put
        all calibrations in the same file, although each part of the response
        function will have its own logical name. We can thus easily modify
        the script to put each calibration information in a separate file.

        Parameters:
         version - Filename version
        Keywords:
         split   - Split IRF over several files?
         clobber - Overwrite existing files?
        """
        # Set calibrate file names
        if split:
            self._ea_file    = "ea_"+version+".fits"
            self._psf_file   = "psf_"+version+".fits"
            self._edisp_file = "edisp_"+version+".fits"
            self._bgd_file   = "bgd_"+version+".fits"
        else:
            self._ea_file    = "irf_"+version+".fits"
            self._psf_file   = "irf_"+version+".fits"
            self._edisp_file = "irf_"+version+".fits"
            self._bgd_file   = "irf_"+version+".fits"

        # Open calibration database index
        self._cif = gammalib.GFits(self._base_path+"/caldb.indx", True)

        # If file has no CIF extension than create it now
        try:
            self._hdu_cif = self._cif.table("CIF")
        except:
            self._create_cif()
            self._hdu_cif = self._cif.table("CIF")

        # Set filenames
        ea_filename    = self._ea_path+"/"+self._ea_file
        psf_filename   = self._psf_path+"/"+self._psf_file
        edisp_filename = self._edisp_path+"/"+self._edisp_file
        bgd_filename   = self._bgd_path+"/"+self._bgd_file

        # Open files
        if split:
            self._ea    = gammalib.GFits(ea_filename, True)
            self._psf   = gammalib.GFits(psf_filename, True)
            self._edisp = gammalib.GFits(edisp_filename, True)
            self._bgd   = gammalib.GFits(bgd_filename, True)
        else:
            self._irf   = gammalib.GFits(ea_filename, True)
            self._ea    = self._irf
            self._psf   = self._irf
            self._edisp = self._irf
            self._bgd   = self._irf

        # Open HDUs
        self._hdu_ea    = self._open_hdu(self._ea, "EFFECTIVE AREA",
                                         self._ea_name, self._ea_doc)
        self._hdu_psf   = self._open_hdu(self._psf, "POINT SPREAD FUNCTION",
                                         self._psf_name, self._psf_doc)
        self._hdu_edisp = self._open_hdu(self._edisp, "ENERGY DISPERSION",
                                         self._edisp_name, self._edisp_doc)
        self._hdu_bgd   = self._open_hdu(self._bgd, "BACKGROUND",
                                         self._bgd_name, self._bgd_doc)

        # Return
        return

    def _close(self):
        """
        Close calibration FITS files.
        """
        # Add information to CIF. We do this now as before we did not
        # necessarily have all the information at hand (in particular
        # about the boundaries)
        self._add_cif_info()

        # Close CIF
        self._close_file(self._cif)

        # If IRF file exists then close it now. All file pointers will
        # be set to None
        if self._irf != None:
            self._close_file(self._irf)
            self._ea    = None
            self._psf   = None
            self._edisp = None
            self._bgd   = None

        # ... otherwise we have split files, so we have to close them all
        else:
            self._close_file(self._ea)
            self._close_file(self._psf)
            self._close_file(self._edisp)
            self._close_file(self._bgd)

        # Return
        return

    def _close_file(self, file):
        """
        Close one file.
        """
        # If file is open then close it
        if file != None:
            file.save(True)
            file.close()
            file = None

        # Return
        return

    def _open_hdu(self, fits, extname, name, doc):
        """
        Open HDU
        
        If HDU does not exist then create one.

        Parameters:
            fits:    FITS file
            extname: Extension name
            name:    Name string
            doc:     Document string

        Returns:
            FITS table
        """
        # Create table if it does not yet exist
        if not fits.contains(extname):

            # Create binary table
            table = gammalib.GFitsBinTable()

            # Set extension name
            table.extname(extname)

            # Set OGIP keywords
            self._set_ogip_keywords(table, doc, ["RESPONSE", name])

            # Append table to FITS file
            fits.append(table)

        # Return FITS table
        return fits.table(extname)

    def _set_ogip_keywords(self, hdu, hdudoc, hduclas):
        """
        Set standard OGIP keywords for extension.

        Parameters:
            hdu:     Header Data Unit.
            hdudoc:  Documentation reference string
            hduclas: Array of HDUCLAS fields
        """
        # Set UTC date of file creation
        utc = datetime.utcnow().isoformat()[:19]

        # Set keywords
        hdu.card("ORIGIN", "IRAP", "Name of organization making this file")
        hdu.card("DATE", utc, "File creation date (YYYY-MM-DDThh:mm:ss UTC)")
        hdu.card("TELESCOP", self._cal_tel, "Name of telescope")
        hdu.card("INSTRUME", self._cal_inst, "Name of instrument")
        hdu.card("DETNAM", self._cal_det, "Name of detector")
        hdu.card("HDUCLASS", "OGIP", "HDU class")
        hdu.card("HDUDOC", hdudoc, "HDU documentation")
        for i, item in enumerate(hduclas):
            key = "HDUCLAS%d" % (i+1)
            hdu.card(key, item, "HDU class")
        hdu.card("HDUVERS", "1.0.0", "HDU version")

        # Return
        return

    def _create_cif(self):
        """
        Create Calibration Index File (CIF) extension in CIF FITS file.

        Parameters:
         None
        Keywords:
         None
        """
        # Create binary table
        table = gammalib.GFitsBinTable()

        # Set boundary

        # Attach columns. Reference: CAL/GEN/92-008
        table.append(gammalib.GFitsTableStringCol("TELESCOP", 0, 10))
        table.append(gammalib.GFitsTableStringCol("INSTRUME", 0, 10))
        table.append(gammalib.GFitsTableStringCol("DETNAM", 0, 20))
        table.append(gammalib.GFitsTableStringCol("FILTER", 0, 10))
        table.append(gammalib.GFitsTableStringCol("CAL_DEV", 0, 20))
        table.append(gammalib.GFitsTableStringCol("CAL_DIR", 0, 70))
        table.append(gammalib.GFitsTableStringCol("CAL_FILE", 0, 40))
        table.append(gammalib.GFitsTableStringCol("CAL_CLAS", 0, 3))
        table.append(gammalib.GFitsTableStringCol("CAL_DTYP", 0, 4))
        table.append(gammalib.GFitsTableStringCol("CAL_CNAM", 0, 20))
        table.append(gammalib.GFitsTableStringCol("CAL_CBD", 0, 70, 9))
        table.append(gammalib.GFitsTableShortCol("CAL_XNO", 0))
        table.append(gammalib.GFitsTableStringCol("CAL_VSD", 0, 10))
        table.append(gammalib.GFitsTableStringCol("CAL_VST", 0, 8))
        table.append(gammalib.GFitsTableDoubleCol("REF_TIME", 0))
        table.append(gammalib.GFitsTableShortCol("CAL_QUAL", 0))
        table.append(gammalib.GFitsTableStringCol("CAL_DATE", 0, 8))
        table.append(gammalib.GFitsTableStringCol("CAL_DESC", 0, 70))

        # Set keywords. Reference: CAL/GEN/92-008
        table.extname("CIF")
        table.card("CIFVERSN", "1992a", "Version of CIF format")

        # Attach table to FITS file. Note that at this moment the
        # FITS table gets connected to the FITS file. Yet since
        # nothing was yet written to the FITS file, we cannot read
        # anything from it.
        self._cif.append(table)

        # Return
        return

    def _add_cif_info(self):
        """
        Add information to CIF extension.
        """
        # Append 4 rows to CIF extension
        self._hdu_cif.append_rows(4)

        # Add generic information for these 4 rows
        for i in range(4):

            # Set row number
            row = i + self._hdu_cif.nrows() - 4

            # Set element
            self._hdu_cif["TELESCOP"][row] = self._cal_tel
            self._hdu_cif["INSTRUME"][row] = self._cal_inst
            self._hdu_cif["DETNAM"][row]   = self._cal_det
            self._hdu_cif["FILTER"][row]   = self._cal_flt
            self._hdu_cif["CAL_DEV"][row]  = "ONLINE"
            self._hdu_cif["CAL_CLAS"][row] = self._cal_class
            self._hdu_cif["CAL_DTYP"][row] = self._cal_type
            self._hdu_cif["CAL_VSD"][row]  = self._val_date
            self._hdu_cif["CAL_VST"][row]  = self._val_time
            self._hdu_cif["REF_TIME"][row] = self._ref_time
            self._hdu_cif["CAL_QUAL"][row] = self._cal_qual
            self._hdu_cif["CAL_DATE"][row] = self._cal_date

        # Add effective area information
        row = self._hdu_cif.nrows() - 4
        self._hdu_cif["CAL_DIR"][row]   = self._ea_dir
        self._hdu_cif["CAL_FILE"][row]  = self._ea_file
        self._hdu_cif["CAL_CNAM"][row]  = self._ea_name
        self._hdu_cif["CAL_DESC"][row]  = self._ea_desc
        self._hdu_cif["CAL_XNO"][row]   = 1
        for n in range(9):
            if n >= len(self._ea_bounds):
                self._hdu_cif["CAL_CBD"][row, n] = "NONE"
            else:
                self._hdu_cif["CAL_CBD"][row, n] = self._ea_bounds[n]

        # Add point spread function information
        row = self._hdu_cif.nrows() - 3
        self._hdu_cif["CAL_DIR"][row]   = self._psf_dir
        self._hdu_cif["CAL_FILE"][row]  = self._psf_file
        self._hdu_cif["CAL_CNAM"][row]  = self._psf_name
        self._hdu_cif["CAL_DESC"][row]  = self._psf_desc
        self._hdu_cif["CAL_XNO"][row]   = 1
        for n in range(9):
            if n >= len(self._psf_bounds):
                self._hdu_cif["CAL_CBD"][row, n] = "NONE"
            else:
                self._hdu_cif["CAL_CBD"][row, n] = self._psf_bounds[n]

        # Add energy dispersion information
        row = self._hdu_cif.nrows() - 2
        self._hdu_cif["CAL_DIR"][row]   = self._edisp_dir
        self._hdu_cif["CAL_FILE"][row]  = self._edisp_file
        self._hdu_cif["CAL_CNAM"][row]  = self._edisp_name
        self._hdu_cif["CAL_DESC"][row]  = self._edisp_desc
        self._hdu_cif["CAL_XNO"][row]   = 1
        for n in range(9):
            if n >= len(self._edisp_bounds):
                self._hdu_cif["CAL_CBD"][row, n] = "NONE"
            else:
                self._hdu_cif["CAL_CBD"][row, n] = self._edisp_bounds[n]

        # Add background information
        row = self._hdu_cif.nrows() - 1
        self._hdu_cif["CAL_DIR"][row]   = self._bgd_dir
        self._hdu_cif["CAL_FILE"][row]  = self._bgd_file
        self._hdu_cif["CAL_CNAM"][row]  = self._bgd_name
        self._hdu_cif["CAL_DESC"][row]  = self._bgd_desc
        self._hdu_cif["CAL_XNO"][row]   = 1
        for n in range(9):
            if n >= len(self._bgd_bounds):
                self._hdu_cif["CAL_CBD"][row, n] = "NONE"
            else:
                self._hdu_cif["CAL_CBD"][row, n] = self._bgd_bounds[n]

        # Return
        return

    def _set_cif_keywords(self, hdu, name, bounds, desc):
        """
        Set standard CIF keywords for extension.
        """
        # Set standard CIF keywords
        hdu.card("CSYSNAME", "XMA_POL", "")
        hdu.card("CCLS0001", self._cal_class, "Calibration class")
        hdu.card("CDTP0001", self._cal_type, "Calibration type")
        hdu.card("CCNM0001", name, "Calibration name")

        # Set boundary keywords
        for n in range(9):
            keyname = "CBD%d0001" % (n+1)
            if n >= len(bounds):
                value = "NONE"
            else:
                value = bounds[n]
            hdu.card(keyname, value, "Calibration boundary")

        # Set validity keywords
        hdu.card("CVSD0001", self._val_date, "Calibration validity start date (UTC)")
        hdu.card("CVST0001", self._val_time, "Calibration validity start time (UTC)")
        hdu.card("CDEC0001", desc, "Calibration description")

        # Return
        return


    # ROOT specific private members
    def _root2caldb(self):
        """
        Translate ROOT to CALDB information.

        Parameters:
         filename - ROOT 2D performance filename.
        """
        # Open ROOT performance file
        file = TFile(self["infile"].filename().url())

        # Create effective area
        self._root2ea(file, rebin=self["rebin"].boolean(),
                            eascale=self["eascale"].real())

        # Create point spread function
        self._root2psf(file, self["psftype"].string())

        # Create energy dispersion
        self._root2edisp(file)

        # Create background
        #self._root2bgd(file)
        self._root2bgd3D(file, bgdscale=self["bgdscale"].real())

        # Return
        return

    def _root2ea(self, file, rebin=False, eascale=1.0):
        """
        Translate ROOT to CALDB effective area extension. The following ROOT
        histograms are used:

        EffectiveAreaEtrue_offaxis -> EFFAREA
        EffectiveArea_offaxis      -> EFFAREA_RECO

        Parameters:
            file: ROOT file.

        Keywords:
            rebin:   Rebin Etrue histogram (useful for Prod1 IFAE runs)
            eascale: Effective area scaling factor
        """
        # Continue only if effective area HDU has been opened
        if self._hdu_ea != None:

            # Allocate ROOT 2D arrays
            etrue = TH2F()
            ereco = TH2F()
            file.GetObject("EffectiveAreaEtrue_offaxis", etrue)
            file.GetObject("EffectiveArea_offaxis",      ereco)

            # Rebin etrue histogram
            if rebin:
                etrue.RebinX(10)
                neng    = etrue.GetXaxis().GetNbins()
                noffset = etrue.GetYaxis().GetNbins()
                for ioff in range(noffset):
                    for ieng in range(neng):
                        value = etrue.GetBinContent(ieng+1,ioff+1) / 10.0
                        etrue.SetBinContent(ieng+1,ioff+1,value)

            # Scale histograms (if needed)
            if eascale != 1.0:
                neng    = etrue.GetXaxis().GetNbins()
                noffset = etrue.GetYaxis().GetNbins()
                for ioff in range(noffset):
                    for ieng in range(neng):
                        value = etrue.GetBinContent(ieng+1,ioff+1) * eascale
                        etrue.SetBinContent(ieng+1,ioff+1,value)
                neng    = ereco.GetXaxis().GetNbins()
                noffset = ereco.GetYaxis().GetNbins()
                for ioff in range(noffset):
                    for ieng in range(neng):
                        value = ereco.GetBinContent(ieng+1,ioff+1) * eascale
                        ereco.SetBinContent(ieng+1,ioff+1,value)

            # Write boundaries (use Ereco boundaries)
            bounds = self._make_2D(ereco, self._hdu_ea, None, "m2")
            for b in bounds:
                self._ea_bounds.append(b)
            self._set_cif_keywords(self._hdu_ea, self._ea_name, \
                                   self._ea_bounds, self._ea_desc)

            # EFFAREA
            self._make_2D(etrue, self._hdu_ea, "EFFAREA", "m2")

            # EFFAREA_RECO
            self._make_2D(ereco, self._hdu_ea, "EFFAREA_RECO", "m2")

        # Return
        return

    def _root2psf(self, file, psftype):
        """
        Translate ROOT to CALDB point spread function extension.
        Parameters:
         file    - ROOT file.
         psftype - PSF type (Gauss, King)
        Keywords:
         None
        """
        # Continue only if point spread function HDU has been opened
        if self._hdu_psf != None:

            # King profile PSF
            if psftype == "King":
                self._root2psf_king(file)
            else:
                self._root2psf_gauss(file)

        # Return
        return

    def _root2psf_gauss(self, file):
        """
        Translate ROOT to CALDB point spread function extension. The following
        ROOT histograms are used:

        1/(2*pi*SIGMA_1) -> SCALE
        AngRes_offaxis -> SIGMA_1 (scaling: 1/0.8)
        0.0 -> AMPL_2
        0.0 -> SIGMA_2
        0.0 -> AMPL_3
        0.0 -> SIGMA_3

        Parameters:
         file - ROOT file.
        Keywords:
         None
        """
        # Continue only if point spread function HDU has been opened
        if self._hdu_psf != None:

            # Allocate ROOT 2D array
            r68  = TH2F()

            # Read 68% containment histogram
            file.GetObject("AngRes_offaxis", r68)
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
            bounds = self._make_2D(r68, self._hdu_psf, None, "deg")
            for b in bounds:
                self._psf_bounds.append(b)
            self._psf_bounds.append("PSF(GAUSS)")
            self._set_cif_keywords(self._hdu_psf, self._psf_name, \
                                   self._psf_bounds, self._psf_desc)

            # SCALE
            self._make_2D(scale, self._hdu_psf, "SCALE", "")

            # SIGMA_1
            self._make_2D(r68, self._hdu_psf, "SIGMA_1", "deg")

            # AMPL_2
            self._make_2D(zero, self._hdu_psf, "AMPL_2", "")

            # SIGMA_2
            self._make_2D(zero, self._hdu_psf, "SIGMA_2", "deg")

            # AMPL_3
            self._make_2D(zero, self._hdu_psf, "AMPL_3", "")

            # SIGMA_3
            self._make_2D(zero, self._hdu_psf, "SIGMA_3", "deg")

        # Return
        return

    def _root2psf_king(self, file):
        """
        Translate ROOT to CALDB point spread function extension. The following
        ROOT histograms are used:

        AngRes_offaxis
        AngRes80_offaxis

        Parameters:
         file - ROOT file.
        Keywords:
         None
        """
        # Continue only if point spread function HDU has been opened
        if self._hdu_psf != None:

            # Allocate ROOT 2D arrays
            r68  = TH2F()
            r80  = TH2F()

            # Read 68% and 80% containment histograms
            file.GetObject("AngRes_offaxis", r68)
            file.GetObject("AngRes80_offaxis", r80)
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
                        iter = 0
                        while True:
                            x     = x1 - f1 * (x1-x2)/(f1-f2)
                            f     = (math.pow(a,x) - 1.0)/(math.pow(b,x) - 1.0) - c
                            iter += 1
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
                            gamma = last_gamma
                            sigma = last_sigma

                        # Store surrent result as last result
                        last_gamma = gamma
                        last_sigma = sigma

                        # Show results on console
                        #print(ioff,ieng,gamma,sigma,x,f,iter)

                    # Set bin contents
                    gamma2D.SetBinContent(ieng+1,ioff+1,gamma)
                    sigma2D.SetBinContent(ieng+1,ioff+1,sigma)

            # Set boundaries
            bounds = self._make_2D(r68, self._hdu_psf, None, "deg")
            for b in bounds:
                self._psf_bounds.append(b)
            self._psf_bounds.append("PSF(KING)")
            self._set_cif_keywords(self._hdu_psf, self._psf_name, \
                                   self._psf_bounds, self._psf_desc)

            # GAMMA
            self._make_2D(gamma2D, self._hdu_psf, "GAMMA", "")

            # SIGMA
            self._make_2D(sigma2D, self._hdu_psf, "SIGMA", "deg")

        # Return
        return

    def _root2edisp(self, file):
        """
        Translate ROOT to CALDB energy dispersion extension. The following ROOT
        histograms are used:

        EestOverEtrue_offaxis  -> MATRIX

        Parameters:
         file - ROOT file.
        Keywords:
         None
        """
        # Continue only if energy dispersion HDU has been opened
        if self._hdu_edisp != None:

            # Allocate ROOT 3D array
            matrix  = TH3F()
            file.GetObject("EestOverEtrue_offaxis", matrix)

            # Set boundaries
            bounds = self._make_3D_migra(matrix, self._hdu_edisp, None, "")
            for b in bounds:
                self._edisp_bounds.append(b)
            self._set_cif_keywords(self._hdu_edisp, self._edisp_name, \
                                   self._edisp_bounds, self._edisp_desc)

            # MATRIX
            self._make_3D_migra(matrix, self._hdu_edisp, "MATRIX", "")

        # Return
        return

    def _root2bgd(self, file):
        """
        Translate ROOT to CALDB background extension. The following ROOT
        histograms are used:

        BGRatePerSqDeg_offaxis -> BGD
        BGRatePerSqDeg_offaxis -> BGD_RECO

        Parameters:
         file - ROOT file.
        Keywords:
         None
        """
        # Continue only if background HDU has been opened
        if self._hdu_bgd != None:

            # Allocate ROOT 2D array
            array = TH2F()
            file.GetObject("BGRatePerSqDeg_offaxis", array)

            # Set boundaries
            bounds = self._make_2D(array, self._hdu_bgd, None, "deg")
            for b in bounds:
                self._bgd_bounds.append(b)
            self._set_cif_keywords(self._hdu_bgd, self._bgd_name, \
                                   self._bgd_bounds, self._bgd_desc)

            # BGD
            self._make_2D(array, self._hdu_bgd, "BGD", "")

            # BGD_RECO
            self._make_2D(array, self._hdu_bgd, "BGD_RECO", "")

        # Return
        return

    def _root2bgd3D(self, file, bgdscale=1.0):
        """
        Translate ROOT to CALDB background extension. The following ROOT
        histograms are used:

        BGRatePerSqDeg_offaxis -> BGD

        Parameters:
         file - ROOT file.
        Keywords:
         scale - Background rate scaling factor
        """
        # Continue only if background HDU has been opened
        if self._hdu_bgd != None:

            # Allocate ROOT 2D array
            array = TH2F()
            file.GetObject("BGRatePerSqDeg_offaxis", array)

            # Set boundaries
            bounds = self._make_3D(array, self._hdu_bgd, None, "deg")
            for b in bounds:
                self._bgd_bounds.append(b)
            self._set_cif_keywords(self._hdu_bgd, self._bgd_name, \
                                   self._bgd_bounds, self._bgd_desc)

            # BGD (reconstructed energy)
            self._make_3D(array, self._hdu_bgd, "BGD", "1/s/MeV/sr", scale=bgdscale)

        # Return
        return

    def _append_column(self, hdu, name, unit, axis, log=False):
        """
        Append column to HDU.
        
        Parameters:
            hdu: HDU
            name: Column name
            unit: Column unit
            axis: ROOT histogram axis
        
        Keywords:
            log:  Logarithmic axis
        """
        # Continue only if columns does not yet exist
        if not hdu.contains(name):

            # Get number of axis bins
            nbins = axis.GetNbins()

            # Append column and set column unit
            hdu.append(gammalib.GFitsTableFloatCol(name, 1, nbins))
            hdu[name].unit(unit)

            # Do we have a LO axis? If not we have a HI axis
            if name[-2:] == "LO":
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

        # Return
        return

    def _make_2D(self, array, hdu, name, unit, scale=1.0):
        """
        Make 2D matrix as function of energy and offset angle from a
        ROOT 2D histogram. If the HDU has already the energy and
        offset angle columns, this method will simply add another data
        column.
        If name==None, the method will not append any data column.

        Parameters:
         array - ROOT 2D histogram.
         hdu   - FITS HDU.
         name  - Data column name.
         unit  - Data unit.
        Keywords:
         scale - Scaling factor for histogram values.
        """
        # Extract energy and offset angle vectors
        energies = array.GetXaxis()
        offsets  = array.GetYaxis()
        neng     = energies.GetNbins()
        noffset  = offsets.GetNbins()

        # Append axis columns to HDU
        self._append_column(hdu, "ENERG_LO", "TeV", energies, log=True)
        self._append_column(hdu, "ENERG_HI", "TeV", energies, log=True)
        self._append_column(hdu, "THETA_LO", "deg", offsets)
        self._append_column(hdu, "THETA_HI", "deg", offsets)

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
        bd_eng = "ENERG(%.4f-%.2f)TeV" % \
                 (pow(10.0, energies.GetBinLowEdge(1)), \
                  pow(10.0, energies.GetBinUpEdge(neng)))
        bd_off = "THETA(%.2f-%.2f)deg" % \
                 (offsets.GetBinLowEdge(1), \
                  offsets.GetBinUpEdge(noffset))
        bd_phi = "PHI(0-360)deg"
        bounds = [bd_eng, bd_off, bd_phi]

        # Return boundary information
        return bounds

    def _make_3D(self, array, hdu, name, unit, scale=1.0):
        """
        Make 3D cube as function of DETX, DETY and energy from a
        ROOT 2D histogram. If the HDU has already the energy and
        offset angle columns, this method will simply add another data
        column.
        If name==None, the method will not append any data column.

        Parameters:
         array - ROOT 2D histogram.
         hdu   - FITS HDU.
         name  - Data column name.
         unit  - Data unit.
        Keywords:
         scale - Scaling factor for histogram values.
        """
        # Set constants
        deg2sr = 0.01745329*0.01745329

        # Extract energy and offset angle vectors
        energies = array.GetXaxis()
        neng     = energies.GetNbins()
        offsets  = array.GetYaxis()
        noffset  = offsets.GetNbins()
        ewidth   = [] # in MeV
        for ieng in range(neng):
            ewidth.append(pow(10.0, energies.GetBinUpEdge(ieng+1) +6.0) - \
                          pow(10.0, energies.GetBinLowEdge(ieng+1)+6.0))

        # Build DETX and DETY axes
        ndet     = array.GetYaxis().GetNbins()
        ndets    = 2*ndet
        det_max  = array.GetYaxis().GetBinUpEdge(ndet)
        det_bin  = det_max/float(ndet)
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

        # Append DETX_LO column
        if not hdu.contains("DETX_LO"):
            hdu.append(gammalib.GFitsTableFloatCol("DETX_LO", 1, ndets))
            hdu["DETX_LO"].unit("deg")
            for i in range(ndets):
                hdu["DETX_LO"][0,i] = dets_lo[i]

        # Append DETX_HI column
        if not hdu.contains("DETX_HI"):
            hdu.append(gammalib.GFitsTableFloatCol("DETX_HI", 1, ndets))
            hdu["DETX_HI"].unit("deg")
            for i in range(ndets):
                hdu["DETX_HI"][0,i] = dets_hi[i]

        # Append DETY_LO column
        if not hdu.contains("DETY_LO"):
            hdu.append(gammalib.GFitsTableFloatCol("DETY_LO", 1, ndets))
            hdu["DETY_LO"].unit("deg")
            for i in range(ndets):
                hdu["DETY_LO"][0,i] = dets_lo[i]

        # Append DETY_HI column
        if not hdu.contains("DETY_HI"):
            hdu.append(gammalib.GFitsTableFloatCol("DETY_HI", 1, ndets))
            hdu["DETY_HI"].unit("deg")
            for i in range(ndets):
                hdu["DETY_HI"][0,i] = dets_hi[i]

        # Append ENERG_LO and ENERG_HI columns
        self._append_column(hdu, "ENERG_LO", "TeV", energies, log=True)
        self._append_column(hdu, "ENERG_HI", "TeV", energies, log=True)

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
                        ioff  = offsets.FindBin(theta)
                        if ioff > 0 and ioff <= noffset:
                            binsq = array.GetBinContent(ieng+1,ioff)
                            if binsq >= 0:
                                value = binsq / deg2sr / (ewidth[ieng])
                                hdu[name][0,index] = value * scale

        # Debugging
        if False:
            for ieng in range(neng):
                sum = 0.0
                for ioff in range(noffset):
                    index = ieng + ioff * neng
                    binsq = array.GetBinContent(ieng+1,ioff+1)
                    value = binsq / deg2sr / ewidth[ieng]
                    sum  += value
                    if ioff == 0:
                        peak = value
                        peaksq = binsq
                print(pow(10.0, energies.GetBinLowEdge(ieng+1))+ \
                      pow(10.0, energies.GetBinUpEdge(ieng+1)), sum, peak, peaksq)

        # Collect boundary information
        bd_detx = "DETX(%.2f-%.2f)deg" % (-det_max, det_max)
        bd_dety = "DETY(%.2f-%.2f)deg" % (-det_max, det_max)
        bd_eng  = "ENERG(%.4f-%.2f)TeV" % \
                  (pow(10.0, energies.GetBinLowEdge(1)),
                   pow(10.0, energies.GetBinUpEdge(neng)))
        bounds  = [bd_detx, bd_dety, bd_eng]

        # Return boundary information
        return bounds

    def _make_3D_migra(self, array, hdu, name, unit, scale=1.0):
        """
        Make 3D cube as function of ETRUE, MIGRA and THETA form ROOT 3D
        histogram. If the HDU has already the columns, this method will
        simply add another data column.
        If name==None, the method will not append any data column.

        Parameters:
         array - ROOT 3D histogram.
         hdu   - FITS HDU.
         name  - Data column name.
         unit  - Data unit.
        Keywords:
         scale - Scaling factor for histogram values.
        """
        # Extract Etrue, Eobs/Etrue and offset angle vectors
        etrue    = array.GetXaxis()
        netrue   = etrue.GetNbins()
        migra    = array.GetYaxis()
        nmigra   = migra.GetNbins()
        offsets  = array.GetZaxis()
        noffset  = offsets.GetNbins()
        ewidth   = [] # in MeV
        for ieng in range(netrue):
            ewidth.append(pow(10.0, etrue.GetBinUpEdge(ieng+1) +6.0) - \
                          pow(10.0, etrue.GetBinLowEdge(ieng+1)+6.0))

        # Append axis columns to HDU
        self._append_column(hdu, "ETRUE_LO", "TeV", etrue, log=True)
        self._append_column(hdu, "ETRUE_HI", "TeV", etrue, log=True)
        self._append_column(hdu, "MIGRA_LO", "", migra)
        self._append_column(hdu, "MIGRA_HI", "", migra)
        self._append_column(hdu, "THETA_LO", "deg", offsets)
        self._append_column(hdu, "THETA_HI", "deg", offsets)

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
        bd_eng   = "ETRUE(%.4f-%.2f)TeV" % \
                   (pow(10.0, etrue.GetBinLowEdge(1)),
                    pow(10.0, etrue.GetBinUpEdge(netrue)))
        bd_migra = "MIGRA(%.3f-%.3f)" % \
                   (migra.GetBinLowEdge(1),
                    migra.GetBinUpEdge(nmigra))
        bd_off   = "THETA(%.2f-%.2f)deg" % \
                   (offsets.GetBinLowEdge(1),
                    offsets.GetBinUpEdge(noffset))
        bd_phi   = "PHI(0-360)deg"
        bounds   = [bd_eng, bd_migra, bd_off, bd_phi]

        # Return boundary information
        return bounds


    # Public methods
    def run(self):
        """
        Run the script.
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

        # Initialise IRF
        self._initialise()

        # Open calibration file
        self._open("file")

        # Translate ROOT to CALDB information
        if _has_root:
            self._root2caldb()

        # Close calibration file
        self._close()

        # Return
        return

    def execute(self):
        """
        Execute the script.
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
