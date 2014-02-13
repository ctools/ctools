#! /usr/bin/env python
# ==========================================================================
# This script generates IRFs in the CALDB format of HEASARC using ROOT 2D
# performance files.
# -------------------------------------------------------------------------
# Copyright (C) 2011-2014 Juergen Knoedlseder
# -------------------------------------------------------------------------
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
# -------------------------------------------------------------------------
# Requirements:
# - gammalib
# - pyroot
# ==========================================================================
from ROOT import TFile, TH1F, TH2F
from gammalib import *
from datetime import datetime
import math
import sys
import glob
import os


# =========== #
# CalDB class #
# =========== #
class caldb():
    """
    """
    def __init__(self, inst, obsid, path=""):
        """
        CALDB constructor. The constructor will create one CALDB entry for
        an observation.

        Parameters:
         inst  - Instrument string.
         obsid - Observation ID string.
        Keywords:
         path  - CALDB Root path.
        """
        # Store root path
        self.path = path
        
        # Initialise some members
        self.cif     = None
        self.irf     = None
        self.ea      = None
        self.psf     = None
        self.edisp   = None
        self.bgd     = None
        self.hdu_cif = None

        # Store UTC date of job execution
        self.utc = datetime.utcnow().isoformat()[:19]
        
        # Set CALDB information
        self.cal_tel      = "CTA"
        self.cal_inst     = inst.upper()
        self.cal_obsid    = obsid
        self.cal_det      = "NONE"
        self.cal_flt      = "NONE"
        self.cal_class    = "BCF"
        self.cal_type     = "DATA"
        self.cal_qual     = 0            # 0=good, 1=bad, 2=dubious, ...
        self.cal_date     = "14/01/30"
        self.val_date     = "2014-01-30"
        self.val_time     = "00:00:00"
        self.ref_time     = 51544.0
        self.cal_name     = "NAME("+obsid+")"
        self.cal_version  = "VERSION(Prod2)"
        self.cal_cut      = "CLASS(BEST)"
        self.cal_analysis = "ANALYSIS(DESY)"
        #
        self.ea_name      = "EFF_AREA"
        self.ea_doc       = "CAL/GEN/92-019"
        self.ea_bounds    = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]
        self.ea_desc      = "CTA effective area"
        #
        self.psf_name     = "RPSF"
        self.psf_doc      = "CAL/GEN/92-020"
        self.psf_bounds   = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]
        self.psf_desc     = "CTA point spread function"
        #
        self.edisp_name   = "EDISP"
        self.edisp_doc    = "???"
        self.edisp_bounds = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]
        self.edisp_desc   = "CTA energy dispersion"
        #
        self.bgd_name     = "BGD"
        self.bgd_doc      = "???"
        self.bgd_bounds   = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]
        self.bgd_desc     = "CTA background"
            
        # Create directory structure
        self.make_dirs()
        
        # Return
        return
    
    def __del__(self):
        """
        Destructor. Close all FITS files.
        """
        # Close calibration
        self.close()

        # Return
        return
    
    def make_dirs(self):
        """
        Generate CALDB directory structure for one observation identifier.
        The structure is given by
        
            data/<tel>/<inst>/bcf/<obsid>
        
        where <tel> is "cta" and <inst> is the instrument specified in the
        CALDB constructor (the instrument may be used for different array
        configurations).
        
        Parameters:
         None
        Keywords:
         None
        """
        # Set base path
        self.base_dir  = "data"
        self.base_dir += "/"+self.cal_tel.lower()
        self.base_dir += "/"+self.cal_inst.lower()

        # Set directory names
        self.bcf_dir   = self.base_dir+"/bcf"
        self.obs_dir   = self.bcf_dir+"/"+self.cal_obsid
        self.ea_dir    = self.obs_dir
        self.psf_dir   = self.obs_dir
        self.edisp_dir = self.obs_dir
        self.bgd_dir   = self.obs_dir
        
        # Optionally prefix path
        if len(self.path) > 0:
            self.base_path  = self.path+"/"+self.base_dir
            self.bcf_path   = self.path+"/"+self.bcf_dir
            self.obs_path   = self.path+"/"+self.obs_dir
            self.ea_path    = self.path+"/"+self.ea_dir
            self.psf_path   = self.path+"/"+self.psf_dir
            self.edisp_path = self.path+"/"+self.edisp_dir
            self.bgd_path   = self.path+"/"+self.bgd_dir
        else:
            self.base_path  = self.base_dir
            self.bcf_path   = self.bcf_dir
            self.obs_path   = self.obs_dir
            self.ea_path    = self.ea_dir
            self.psf_path   = self.psf_dir
            self.edisp_path = self.edisp_dir
            self.bgd_path   = self.bgd_dir

        # Create directory structure
        if not os.path.isdir(self.ea_path):
            os.makedirs(self.ea_path)
        if not os.path.isdir(self.psf_path):
            os.makedirs(self.psf_path)
        if not os.path.isdir(self.edisp_path):
            os.makedirs(self.edisp_path)
        if not os.path.isdir(self.bgd_path):
            os.makedirs(self.bgd_path)
    
        # Return
        return
    
    def open(self, version, split=False, clobber=True):
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
            self.ea_file    = "ea_"+version+".fits"
            self.psf_file   = "psf_"+version+".fits"
            self.edisp_file = "edisp_"+version+".fits"
            self.bgd_file   = "bgd_"+version+".fits"
        else:
            self.ea_file    = "irf_"+version+".fits"
            self.psf_file   = "irf_"+version+".fits"
            self.edisp_file = "irf_"+version+".fits"
            self.bgd_file   = "irf_"+version+".fits"
        
        # Open calibration database index
        self.cif = GFits(self.base_path+"/caldb.indx", True)
        
        # If file has no CIF extension than create it now
        try:
            self.hdu_cif = self.cif.table("CIF")
        except:
            self.create_cif()
            self.hdu_cif = self.cif.table("CIF")
        
        # Set filenames
        ea_filename    = self.ea_path+"/"+self.ea_file
        psf_filename   = self.psf_path+"/"+self.psf_file
        edisp_filename = self.edisp_path+"/"+self.edisp_file
        bgd_filename   = self.bgd_path+"/"+self.bgd_file

        # Open files
        if split:
            self.ea    = GFits(ea_filename, True)
            self.psf   = GFits(psf_filename, True)
            self.edisp = GFits(edisp_filename, True)
            self.bgd   = GFits(bgd_filename, True)
        else:
            self.irf   = GFits(ea_filename, True)
            self.ea    = self.irf
            self.psf   = self.irf
            self.edisp = self.irf
            self.bgd   = self.irf
        
        # Open HDUs
        self.open_ea()
        self.open_psf()
        self.open_edisp()
        self.open_bgd()
        
        # Return
        return

    def close(self):
        """
        Close calibration FITS files.

        Parameters:
         None
        Keywords:
         None
        """
        # Add information to CIF. We do this now as before we did not
        # necessarily have all the information at hand (in particular
        # about the boundaries)
        self.add_cif_info()
        
        # Close CIF
        if self.cif != None:
            self.cif.save(True)
            self.cif.close()
            self.cif = None

        # If IRF file exists then close it now. All file pointers will
        # be set to None
        if self.irf != None:
            self.irf.save(True)
            self.irf.close()
            self.irf   = None
            self.ea    = None
            self.psf   = None
            self.edisp = None
            self.bgd   = None
        
        # ... otherwise we have split files, so we have to close them
        # all
        else:
        
            # Close effective area file
            if self.ea != None:
                self.ea.save(True)
                self.ea.close()
                self.ea = None

            # Close point spread function file
            if self.psf != None:
                self.psf.save(True)
                self.psf.close()
                self.psf = None

            # Close energy dispersion file
            if self.edisp != None:
                self.edisp.save(True)
                self.edisp.close()
                self.edisp = None

            # Close background file
            if self.bgd != None:
                self.bgd.save(True)
                self.bgd.close()
                self.bgd = None
        
        # Return
        return
    
    def create_cif(self):
        """
        Create Calibration Index File (CIF) extension in CIF FITS file.

        Parameters:
         None
        Keywords:
         None
        """
        # Create binary table
        table = GFitsBinTable()
        
        # Set boundary
        
        # Attach columns. Reference: CAL/GEN/92-008
        table.append(GFitsTableStringCol("TELESCOP", 0, 10))
        table.append(GFitsTableStringCol("INSTRUME", 0, 10))
        table.append(GFitsTableStringCol("DETNAM", 0, 20))
        table.append(GFitsTableStringCol("FILTER", 0, 10))
        table.append(GFitsTableStringCol("CAL_DEV", 0, 20))
        table.append(GFitsTableStringCol("CAL_DIR", 0, 70))
        table.append(GFitsTableStringCol("CAL_FILE", 0, 40))
        table.append(GFitsTableStringCol("CAL_CLAS", 0, 3))
        table.append(GFitsTableStringCol("CAL_DTYP", 0, 4))
        table.append(GFitsTableStringCol("CAL_CNAM", 0, 20))
        table.append(GFitsTableStringCol("CAL_CBD", 0, 70, 9))
        table.append(GFitsTableShortCol("CAL_XNO", 0))
        table.append(GFitsTableStringCol("CAL_VSD", 0, 10))
        table.append(GFitsTableStringCol("CAL_VST", 0, 8))
        table.append(GFitsTableDoubleCol("REF_TIME", 0))
        table.append(GFitsTableShortCol("CAL_QUAL", 0))
        table.append(GFitsTableStringCol("CAL_DATE", 0, 8))
        table.append(GFitsTableStringCol("CAL_DESC", 0, 70))
        
        # Set keywords. Reference: CAL/GEN/92-008
        table.extname("CIF")
        table.card("CIFVERSN", "1992a", "Version of CIF format")
        
        # Attach table to FITS file. Note that at this moment the
        # FITS table gets connected to the FITS file. Yet since
        # nothing was yet written to the FITS file, we cannot read
        # anything from it.
        self.cif.append(table)
        
        # Return
        return
    
    def add_cif_info(self):
        """
        Add information to CIF extension.

        Parameters:
         None
        Keywords:
         None
        """
        # Append 4 rows to CIF extension
        self.hdu_cif.append_rows(4)
        
        # Add generic information for these 4 rows
        for i in range(4):
        
            # Set row number
            row = i+self.hdu_cif.nrows()-4
            
            # Set element
            self.hdu_cif["TELESCOP"][row] = self.cal_tel
            self.hdu_cif["INSTRUME"][row] = self.cal_inst
            self.hdu_cif["DETNAM"][row]   = self.cal_det 
            self.hdu_cif["FILTER"][row]   = self.cal_flt
            self.hdu_cif["CAL_DEV"][row]  = "ONLINE"
            self.hdu_cif["CAL_CLAS"][row] = self.cal_class
            self.hdu_cif["CAL_DTYP"][row] = self.cal_type
            self.hdu_cif["CAL_VSD"][row]  = self.val_date
            self.hdu_cif["CAL_VST"][row]  = self.val_time
            self.hdu_cif["REF_TIME"][row] = self.ref_time
            self.hdu_cif["CAL_QUAL"][row] = self.cal_qual
            self.hdu_cif["CAL_DATE"][row] = self.cal_date
        
        # Add effective area information
        row = self.hdu_cif.nrows()-4
        self.hdu_cif["CAL_DIR"][row]   = self.ea_dir
        self.hdu_cif["CAL_FILE"][row]  = self.ea_file
        self.hdu_cif["CAL_CNAM"][row]  = self.ea_name
        self.hdu_cif["CAL_DESC"][row]  = self.ea_desc
        self.hdu_cif["CAL_XNO"][row]   = 1
        for n in range(9):
            if n >= len(self.ea_bounds):
                self.hdu_cif["CAL_CBD"][row, n] = "NONE"
            else:
                self.hdu_cif["CAL_CBD"][row, n] = self.ea_bounds[n]

        # Add point spread function information
        row = self.hdu_cif.nrows()-3
        self.hdu_cif["CAL_DIR"][row]   = self.psf_dir
        self.hdu_cif["CAL_FILE"][row]  = self.psf_file
        self.hdu_cif["CAL_CNAM"][row]  = self.psf_name
        self.hdu_cif["CAL_DESC"][row]  = self.psf_desc
        self.hdu_cif["CAL_XNO"][row]   = 1
        for n in range(9):
            if n >= len(self.psf_bounds):
                self.hdu_cif["CAL_CBD"][row, n] = "NONE"
            else:
                self.hdu_cif["CAL_CBD"][row, n] = self.psf_bounds[n]

        # Add energy dispersion information
        row = self.hdu_cif.nrows()-2
        self.hdu_cif["CAL_DIR"][row]   = self.edisp_dir
        self.hdu_cif["CAL_FILE"][row]  = self.edisp_file
        self.hdu_cif["CAL_CNAM"][row]  = self.edisp_name
        self.hdu_cif["CAL_DESC"][row]  = self.edisp_desc
        self.hdu_cif["CAL_XNO"][row]   = 1
        for n in range(9):
            if n >= len(self.edisp_bounds):
                self.hdu_cif["CAL_CBD"][row, n] = "NONE"
            else:
                self.hdu_cif["CAL_CBD"][row, n] = self.edisp_bounds[n]

        # Add background information
        row = self.hdu_cif.nrows()-1
        self.hdu_cif["CAL_DIR"][row]   = self.bgd_dir
        self.hdu_cif["CAL_FILE"][row]  = self.bgd_file
        self.hdu_cif["CAL_CNAM"][row]  = self.bgd_name
        self.hdu_cif["CAL_DESC"][row]  = self.bgd_desc
        self.hdu_cif["CAL_XNO"][row]   = 1
        for n in range(9):
            if n >= len(self.bgd_bounds):
                self.hdu_cif["CAL_CBD"][row, n] = "NONE"
            else:
                self.hdu_cif["CAL_CBD"][row, n] = self.bgd_bounds[n]
        
        # Return
        return

    def open_ea(self):
        """
        Open effective area extension. If no extension was found
        then create one. We do not yet set any data as we don't
        know the dimensions of the table yet.

        Parameters:
         None
        Keywords:
         None
        """
        # Get extension. If it does not exist then create it now.
        try:
            self.hdu_ea = self.ea.table("EFFECTIVE AREA")
        except:
            # Create binary table
            table = GFitsBinTable()
            
            # Set extension name
            table.extname("EFFECTIVE AREA")
            
            # Set OGIP keywords
            self.set_ogip_keywords(table, self.ea_doc, \
                                   ["RESPONSE", self.ea_name])
            
            # Append table
            self.ea.append(table)
            
            # Get extension
            self.hdu_ea = self.ea.table("EFFECTIVE AREA")
        
        # Return
        return

    def open_psf(self):
        """
        Open point spread function extension. If no extension was
        found then create one. We do not yet set any data as we don't
        know the dimensions of the table yet.

        Parameters:
         None
        Keywords:
         None
        """
        # Get extension. If it does not exist then create it now.
        try:
            self.hdu_psf = self.psf.table("POINT SPREAD FUNCTION")
        except:
            # Create binary table
            table = GFitsBinTable()
            
            # Set extension name
            table.extname("POINT SPREAD FUNCTION")
            
            # Set OGIP keywords
            self.set_ogip_keywords(table, self.psf_doc, \
                                   ["RESPONSE", self.psf_name])
            
            # Append table
            self.psf.append(table)

            # Get extension
            self.hdu_psf = self.psf.table("POINT SPREAD FUNCTION")
        
        # Return
        return

    def open_edisp(self):
        """
        Open energy dispersion extension. If no extension was found
        then create one. We do not yet set any data as we don't
        know the dimensions of the table yet.

        Parameters:
         None
        Keywords:
         None
        """
        # Get extension. If it does not exist then create it now.
        try:
            self.hdu_edisp = self.edisp.table("ENERGY DISPERSION")
        except:
            # Create binary table
            table = GFitsBinTable()
            
            # Set extension name
            table.extname("ENERGY DISPERSION")
            
            # Set OGIP keywords
            self.set_ogip_keywords(table, self.edisp_doc, \
                                   ["RESPONSE", self.edisp_name])

            # Append table
            self.edisp.append(table)

            # Get extension
            self.hdu_edisp = self.edisp.table("ENERGY DISPERSION")
        
        # Return
        return

    def open_bgd(self):
        """
        Open background extension. If no extension was found
        then create one. We do not yet set any data as we don't
        know the dimensions of the table yet.

        Parameters:
         None
        Keywords:
         None
        """
        # Get extension. If it does not exist then create it now.
        try:
            self.hdu_bgd = self.bgd.table("BACKGROUND")
        except:
            # Create binary table
            table = GFitsBinTable()
            
            # Set extension name
            table.extname("BACKGROUND")
            
            # Set OGIP keywords
            self.set_ogip_keywords(table, self.bgd_doc, \
                                   ["RESPONSE", self.bgd_name])

            # Append table
            self.bgd.append(table)

            # Get extension
            self.hdu_bgd = self.bgd.table("BACKGROUND")
        
        # Return
        return

    def set_ogip_keywords(self, hdu, hdudoc, hduclas):
        """
        Set standard OGIP keywords for extension.

        Parameters:
         hdu     - Header Data Unit.
         hdudoc  - Documentation reference string
         hduclas - Array of HDUCLAS fields
        Keywords:
         None
        """
        # Set keywords
        hdu.card("ORIGIN", "IRAP", "Name of organization making this file")
        hdu.card("DATE", self.utc, "File creation date (YYYY-MM-DDThh:mm:ss UTC)")
        hdu.card("TELESCOP", self.cal_tel, "Name of telescope")
        hdu.card("INSTRUME", self.cal_inst, "Name of instrument")
        hdu.card("DETNAM", self.cal_det, "Name of detector")
        hdu.card("HDUCLASS", "OGIP", "")
        hdu.card("HDUDOC", hdudoc, "")
        for i, item in enumerate(hduclas):
            key = "HDUCLAS%d" % (i+1)
            hdu.card(key, item, "")
        hdu.card("HDUVERS", "1.0.0", "")
        
        # Return
        return

    def make_2D(self, array, hdu, name, unit, scale=1.0):
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

        # Attach columns to HDU. We do this only if they do not yet
        # exist. This allows adding columns to the same HDU in successive
        # calls of this method
        #
        # ENERG_LO
        if not hdu.contains("ENERG_LO"):
            hdu.append(GFitsTableFloatCol("ENERG_LO", 1, neng))
            hdu["ENERG_LO"].unit("TeV")
            for ieng in range(neng):
                e_lo = pow(10.0, energies.GetBinLowEdge(ieng+1))
                hdu["ENERG_LO"][0,ieng] = e_lo
        #
        # ENERG_HI
        if not hdu.contains("ENERG_HI"):
            hdu.append(GFitsTableFloatCol("ENERG_HI", 1, neng))
            hdu["ENERG_HI"].unit("TeV")
            for ieng in range(neng):
                e_hi = pow(10.0, energies.GetBinUpEdge(ieng+1))
                hdu["ENERG_HI"][0,ieng] = e_hi
        #
        # THETA_LO
        if not hdu.contains("THETA_LO"):
            hdu.append(GFitsTableFloatCol("THETA_LO", 1, noffset))
            hdu["THETA_LO"].unit("deg")
            for ioff in range(noffset):
                o_lo = offsets.GetBinLowEdge(ioff+1)
                hdu["THETA_LO"][0,ioff] = o_lo
        #
        # THETA_LO
        if not hdu.contains("THETA_HI"):
            hdu.append(GFitsTableFloatCol("THETA_HI", 1, noffset))
            hdu["THETA_HI"].unit("deg")
            for ioff in range(noffset):
                o_hi = offsets.GetBinUpEdge(ioff+1)
                hdu["THETA_HI"][0,ioff] = o_hi
        #
        # "NAME"
        if name != None and not hdu.contains(name):
            hdu.append(GFitsTableFloatCol(name, 1, neng*noffset))
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

    def make_3D(self, array, hdu, name, unit, scale=1.0):
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

        # Attach columns to HDU. We do this only if they do not yet
        # exist. This allows adding columns to the same HDU in successive
        # calls of this method
        #
        # DETX_LO
        if not hdu.contains("DETX_LO"):
            hdu.append(GFitsTableFloatCol("DETX_LO", 1, ndets))
            hdu["DETX_LO"].unit("deg")
            for i in range(ndets):
                hdu["DETX_LO"][0,i] = dets_lo[i]
        #
        # DETX_HI
        if not hdu.contains("DETX_HI"):
            hdu.append(GFitsTableFloatCol("DETX_HI", 1, ndets))
            hdu["DETX_HI"].unit("deg")
            for i in range(ndets):
                hdu["DETX_HI"][0,i] = dets_hi[i]
        #
        # DETY_LO
        if not hdu.contains("DETY_LO"):
            hdu.append(GFitsTableFloatCol("DETY_LO", 1, ndets))
            hdu["DETY_LO"].unit("deg")
            for i in range(ndets):
                hdu["DETY_LO"][0,i] = dets_lo[i]
        #
        # DETY_HI
        if not hdu.contains("DETY_HI"):
            hdu.append(GFitsTableFloatCol("DETY_HI", 1, ndets))
            hdu["DETY_HI"].unit("deg")
            for i in range(ndets):
                hdu["DETY_HI"][0,i] = dets_hi[i]
        #
        # ENERG_LO
        if not hdu.contains("ENERG_LO"):
            hdu.append(GFitsTableFloatCol("ENERG_LO", 1, neng))
            hdu["ENERG_LO"].unit("TeV")
            for ieng in range(neng):
                e_lo = pow(10.0, energies.GetBinLowEdge(ieng+1))
                hdu["ENERG_LO"][0,ieng] = e_lo
        #
        # ENERG_HI
        if not hdu.contains("ENERG_HI"):
            hdu.append(GFitsTableFloatCol("ENERG_HI", 1, neng))
            hdu["ENERG_HI"].unit("TeV")
            for ieng in range(neng):
                e_hi = pow(10.0, energies.GetBinUpEdge(ieng+1))
                hdu["ENERG_HI"][0,ieng] = e_hi
        #
        # "NAME"
        if name != None and not hdu.contains(name):
            hdu.append(GFitsTableFloatCol(name, 1, ndets*ndets*neng))
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
        #
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
                  (pow(10.0, energies.GetBinLowEdge(1)), \
                   pow(10.0, energies.GetBinUpEdge(neng)))
        bounds  = [bd_detx, bd_dety, bd_eng]

        # Return boundary information
        return bounds
    
    def root2caldb(self, filename, rebin=False, psftype="Gauss", scale=1.0):
        """
        Translate ROOT to CALDB information.
        
        Parameters:
         filename - ROOT 2D performance filename.
        Keywords:
         rebin - Rebin Etrue histogram (useful for Prod1 IFAE runs)
        """
        # Open ROOT performance file
        file = TFile(filename)

        # Create effective area
        self.root2ea(file, rebin=rebin)
        
        # Create point spread function
        self.root2psf(file, psftype)

        # Create energy dispersion
        self.root2edisp(file)

        # Create background
        #self.root2bgd(file)
        self.root2bgd3D(file, scale=scale)
        
        # Return
        return

    def root2ea(self, file, rebin=False):
        """
        Translate ROOT to CALDB effective area extension. The following ROOT
        histograms are used:

        EffectiveAreaEtrue_offaxis -> EFFAREA
        EffectiveArea_offaxis      -> EFFAREA_RECO

        Parameters:
         file - ROOT file.
        Keywords:
         rebin - Rebin Etrue histogram (useful for Prod1 IFAE runs)
        """
        # Continue only if effective area HDU has been opened
        if self.hdu_ea != None:

            # Allocate ROOT 2D arrays
            etrue = TH2F()
            ereco = TH2F()
            file.GetObject("EffectiveAreaEtrue_offaxis", etrue)
            file.GetObject("EffectiveArea_offaxis",      ereco)
            #print(etrue.GetXaxis().GetBinLowEdge(1), etrue.GetXaxis().GetBinUpEdge(etrue.GetXaxis().GetNbins()))
            #print(ereco.GetXaxis().GetBinLowEdge(1), ereco.GetXaxis().GetBinUpEdge(ereco.GetXaxis().GetNbins()))

            # Rebin etrue histogram
            if rebin:
                etrue.RebinX(10)
                neng    = etrue.GetXaxis().GetNbins()
                noffset = etrue.GetYaxis().GetNbins()
                for ioff in range(noffset):
                    for ieng in range(neng):
                        value = etrue.GetBinContent(ieng+1,ioff+1) / 10.0
                        etrue.SetBinContent(ieng+1,ioff+1,value)

            # Write boundaries (use Ereco boundaries)
            bounds = self.make_2D(ereco, self.hdu_ea, None, "m2")
            for b in bounds:
                self.ea_bounds.append(b)
            self.set_cif_keywords(self.hdu_ea, self.ea_name, \
                                  self.ea_bounds, self.ea_desc)

            # EFFAREA
            self.make_2D(etrue, self.hdu_ea, "EFFAREA", "m2")

            # EFFAREA_RECO
            self.make_2D(ereco, self.hdu_ea, "EFFAREA_RECO", "m2")

        # Return
        return

    def root2psf(self, file, psftype):
        """
        Translate ROOT to CALDB point spread function extension.
        Parameters:
         file    - ROOT file.
         psftype - PSF type (Gauss, King)
        Keywords:
         None
        """
        # Continue only if point spread function HDU has been opened
        if self.hdu_psf != None:

            # King profile PSF
            if psftype == "King":
                self.root2psf_king(file)
            else:
                self.root2psf_gauss(file)

        # Return
        return

    def root2psf_gauss(self, file):
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
        if self.hdu_psf != None:

            # Allocate ROOT 2D array
            r68  = TH2F()

            # Read 68% containment histogram
            file.GetObject("AngRes_offaxis", r68)
            neng    = r68.GetXaxis().GetNbins()
            noffset = r68.GetYaxis().GetNbins()
            
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
            bounds = self.make_2D(r68, self.hdu_psf, None, "deg")
            for b in bounds:
                self.psf_bounds.append(b)
            self.psf_bounds.append("PSF(GAUSS)")
            self.set_cif_keywords(self.hdu_psf, self.psf_name, \
                                  self.psf_bounds, self.psf_desc)

            # SCALE
            self.make_2D(scale, self.hdu_psf, "SCALE", "")

            # SIGMA_1
            self.make_2D(r68, self.hdu_psf, "SIGMA_1", "deg")

            # AMPL_2
            self.make_2D(zero, self.hdu_psf, "AMPL_2", "")

            # SIGMA_2
            self.make_2D(zero, self.hdu_psf, "SIGMA_2", "deg")

            # AMPL_3
            self.make_2D(zero, self.hdu_psf, "AMPL_3", "")

            # SIGMA_3
            self.make_2D(zero, self.hdu_psf, "SIGMA_3", "deg")

        # Return
        return

    def root2psf_king(self, file):
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
        if self.hdu_psf != None:

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
            bounds = self.make_2D(r68, self.hdu_psf, None, "deg")
            for b in bounds:
                self.psf_bounds.append(b)
            self.psf_bounds.append("PSF(KING)")
            self.set_cif_keywords(self.hdu_psf, self.psf_name, \
                                  self.psf_bounds, self.psf_desc)

            # GAMMA
            self.make_2D(gamma2D, self.hdu_psf, "GAMMA", "")

            # SIGMA
            self.make_2D(sigma2D, self.hdu_psf, "SIGMA", "deg")

        # Return
        return

    def root2edisp(self, file):
        """
        Translate ROOT to CALDB energy dispersion extension. The following ROOT
        histograms are used:

        ERes_offaxis  -> ERES
        Ebias_offaxis -> EBIAS

        Parameters:
         file - ROOT file.
        Keywords:
         None
        """
        # Continue only if energy dispersion HDU has been opened
        if self.hdu_edisp != None:

            # Allocate ROOT 2D array
            eres  = TH2F()
            ebias = TH2F()
            file.GetObject("ERes_offaxis", eres)
            file.GetObject("Ebias_offaxis", ebias)

            # Set boundaries
            bounds = self.make_2D(eres, self.hdu_edisp, None, "deg")
            for b in bounds:
                self.edisp_bounds.append(b)
            self.set_cif_keywords(self.hdu_edisp, self.edisp_name, \
                                  self.edisp_bounds, self.edisp_desc)

            # ERES
            self.make_2D(eres, self.hdu_edisp, "ERES", "Ereco/Etrue")

            # EBIAS
            self.make_2D(ebias, self.hdu_edisp, "EBIAS", "Ereco/Etrue")

        # Return
        return

    def root2bgd(self, file):
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
        if self.hdu_bgd != None:

            # Allocate ROOT 2D array
            array = TH2F()
            file.GetObject("BGRatePerSqDeg_offaxis", array)

            # Set boundaries
            bounds = self.make_2D(array, self.hdu_bgd, None, "deg")
            for b in bounds:
                self.bgd_bounds.append(b)
            self.set_cif_keywords(self.hdu_bgd, self.bgd_name, \
                                  self.bgd_bounds, self.bgd_desc)

            # BGD
            self.make_2D(array, self.hdu_bgd, "BGD", "")

            # BGD_RECO
            self.make_2D(array, self.hdu_bgd, "BGD_RECO", "")
            
        # Return
        return

    def root2bgd3D(self, file, scale=1.0):
        """
        Translate ROOT to CALDB background extension. The following ROOT
        histograms are used:

        BGRatePerSqDeg_offaxis -> BGD
        BGRatePerSqDeg_offaxis -> BGD_RECO

        Parameters:
         file - ROOT file.
        Keywords:
         scale - Background rate scaling factor
        """
        # Continue only if background HDU has been opened
        if self.hdu_bgd != None:

            # Allocate ROOT 2D array
            array = TH2F()
            file.GetObject("BGRatePerSqDeg_offaxis", array)

            # Set boundaries
            bounds = self.make_3D(array, self.hdu_bgd, None, "deg")
            for b in bounds:
                self.bgd_bounds.append(b)
            self.set_cif_keywords(self.hdu_bgd, self.bgd_name, \
                                  self.bgd_bounds, self.bgd_desc)

            # BGD
            self.make_3D(array, self.hdu_bgd, "BGD", "1/s/MeV/sr", scale=scale)

            # BGD_RECO
            self.make_3D(array, self.hdu_bgd, "BGD_RECO", "1/s/MeV/sr", scale=scale)
            
        # Return
        return

    def set_cif_keywords(self, hdu, name, bounds, desc):
        """
        Set standard CIF keywords for extension.
        """
        # Set standard CIF keywords
        hdu.card("CSYSNAME", "XMA_POL", "")
        hdu.card("CCLS0001", self.cal_class, "Calibration class")
        hdu.card("CDTP0001", self.cal_type, "Calibration type")
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
        hdu.card("CVSD0001", self.val_date, "Calibration validity start date (UTC)")
        hdu.card("CVST0001", self.val_time, "Calibration validity start time (UTC)")
        hdu.card("CDEC0001", desc, "Calibration description")
        
        # Return
        return


#===================#
# Set test database #
#===================#
def set_test():
    """
    Set Prod1 IFAE database.
    """
    # Set database attributes
    path    = "/project-data/cta/performance/prod1/IFAEOffaxisPerformanceBEI_May2012"
    rebin   = True
    scale   = 1.0

    # Set database content
    db = [{'inst': "e", 'id':   "IFAE20120510_50h",
           'path': path, 'rebin': rebin, 'psftype': "Gauss", 'scale': scale,
           'file': "SubarrayE_IFAE_50hours_20120510_offaxis.root"},
          {'inst': "e", 'id':   "IFAE20120510_50h_King",
           'path': path, 'rebin': rebin, 'psftype': "King", 'scale': scale,
           'file': "SubarrayE_IFAE_50hours_20120510_offaxis.root"}]

    # Return database
    return db


#=========================#
# Set Prod1 IFAE database #
#=========================#
def set_prod1_ifae():
    """
    Set Prod1 IFAE database.
    """
    # Set database attributes
    path    = "/project-data/cta/performance/prod1/IFAEOffaxisPerformanceBEI_May2012"
    rebin   = True
    psftype = "Gauss"
    scale   = 1.0

    # Set database content
    db = [{'inst': "b", 'id':   "IFAE20120510_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "SubarrayB_IFAE_50hours_20120510_offaxis.root"},
          {'inst': "e", 'id':   "IFAE20120510_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "SubarrayE_IFAE_50hours_20120510_offaxis.root"},
          {'inst': "i", 'id':   "IFAE20120510_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "SubarrayI_IFAE_50hours_20120510_offaxis.root"}]

    # Return database
    return db


#=========================#
# Set Prod2 DESY database #
#=========================#
def set_prod2_desy_ignore():
    """
    Set Prod2 DESY database (ignore these files).
    """
    # Set database attributes
    path    = "/project-data/cta/performance/prod2/Performance_DESY_20140128"
    rebin   = False
    psftype = "Gauss"
    scale   = 1.0/180000.0

    # Set database content
    db = [{'inst': "aar", 'id':   "DESY20140105_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0NIM2.prod2-Aar-NS.S.2a.180000s.root"},
          {'inst': "aar", 'id':   "DESY20140105_50h_0deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_0degNIM2.prod2-Aar-NS.S.2a.180000s.root"},
          {'inst': "aar", 'id':   "DESY20140105_50h_180deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_180degNIM2.prod2-Aar-NS.S.2a.180000s.root"},

          {'inst': "aar500", 'id':   "DESY20140105_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0NIM2.prod2-Aar-500m-NS.S.2a.180000s.root"},
          {'inst': "aar500", 'id':   "DESY20140105_50h_0deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_0degNIM2.prod2-Aar-500m-NS.S.2a.180000s.root"},
          {'inst': "aar500", 'id':   "DESY20140105_50h_180deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_180degNIM2.prod2-Aar-500m-NS.S.2a.180000s.root"},

          {'inst': "leoncito", 'id':   "DESY20140105_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0NIM2.prod2-LeoncitoPP-NS.S.2a.180000s.root"},
          {'inst': "leoncito", 'id':   "DESY20140105_50h_0deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_0degNIM2.prod2-LeoncitoPP-NS.S.2a.180000s.root"},
          {'inst': "leoncito", 'id':   "DESY20140105_50h_180deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_180degNIM2.prod2-LeoncitoPP-NS.S.2a.180000s.root"},

          {'inst': "sac", 'id':   "DESY20140105_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0NIM2.prod2-SAC100-NS.S.2a.180000s.root"},
          {'inst': "sac", 'id':   "DESY20140105_50h_0deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_0degNIM2.prod2-SAC100-NS.S.2a.180000s.root"},
          {'inst': "sac", 'id':   "DESY20140105_50h_180deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_180degNIM2.prod2-SAC100-NS.S.2a.180000s.root"},

          {'inst': "spm", 'id':   "DESY20140105_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0NIM2.prod2-SPM-NS.N.2NN.180000s.root"},
          {'inst': "spm", 'id':   "DESY20140105_50h_0deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_0degNIM2.prod2-SPM-NS.N.2NN.180000s.root"},
          {'inst': "spm", 'id':   "DESY20140105_50h_180deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_180degNIM2.prod2-SPM-NS.N.2NN.180000s.root"},

          {'inst': "tenerife", 'id':   "DESY20140105_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0NIM2.prod2-Tenerife-NS.N.2NN.180000s.root"},
          {'inst': "tenerife", 'id':   "DESY20140105_50h_0deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_0degNIM2.prod2-Tenerife-NS.N.2NN.180000s.root"},
          {'inst': "tenerife", 'id':   "DESY20140105_50h_180deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_180degNIM2.prod2-Tenerife-NS.N.2NN.180000s.root"},

          {'inst': "us", 'id':   "DESY20140105_50h",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0NIM2.prod2-US-NS.N.2NN.180000s.root"},
          {'inst': "us", 'id':   "DESY20140105_50h_0deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_0degNIM2.prod2-US-NS.N.2NN.180000s.root"},
          {'inst': "us", 'id':   "DESY20140105_50h_180deg",
           'path': path, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V2.ID0_180degNIM2.prod2-US-NS.N.2NN.180000s.root"}]

    # Return database
    return db


#=========================#
# Set Prod2 DESY database #
#=========================#
def set_prod2_desy():
    """
    Set Prod2 DESY database.
    """
    # Set database attributes
    path_aar = "/project-data/cta/performance/prod2/Performance_DESY_20140131_Aar"
    path_ten = "/project-data/cta/performance/prod2/Performance_DESY_20140131_Tenerife_50h"
    rebin    = False
    psftype  = "Gauss"
    scale    = 1.0

    # Set database content
    db = [{'inst': "aar", 'id':   "DESY20140105_50h",
           'path': path_aar, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0NIM2.prod2-Aar-NS.S.2a.180000s.root"},
          {'inst': "aar", 'id':   "DESY20140105_50h_0deg",
           'path': path_aar, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0_0degNIM2.prod2-Aar-NS.S.2a.180000s.root"},
          {'inst': "aar", 'id':   "DESY20140105_50h_180deg",
           'path': path_aar, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0_180degNIM2.prod2-Aar-NS.S.2a.180000s.root"},

          {'inst': "aar500", 'id':   "DESY20140105_50h",
           'path': path_aar, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0NIM2.prod2-Aar-500m-NS.S.2a.180000s.root"},
          {'inst': "aar500", 'id':   "DESY20140105_50h_0deg",
           'path': path_aar, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0_0degNIM2.prod2-Aar-NS.S.2a.180000s.root"},
          {'inst': "aar500", 'id':   "DESY20140105_50h_180deg",
           'path': path_aar, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0_180degNIM2.prod2-Aar-NS.S.2a.180000s.root"},

          {'inst': "tenerife", 'id':   "DESY20140105_50h",
           'path': path_ten, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0NIM2.prod2-Tenerife-NS.N.2NN.180000s.root"},
          {'inst': "tenerife", 'id':   "DESY20140105_50h_0deg",
           'path': path_ten, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0_0degNIM2.prod2-Tenerife-NS.N.2NN.180000s.root"},
          {'inst': "tenerife", 'id':   "DESY20140105_50h_180deg",
           'path': path_ten, 'rebin': rebin, 'psftype': psftype, 'scale': scale,
           'file': "DESY.d20140105.Erec1.V3.ID0_180degNIM2.prod2-Tenerife-NS.N.2NN.180000s.root"}
         ]

    # Return database
    return db


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    This script translates a CTA ROOT 2D performance files into a CALDB
    compliant CTA response. The script will create the response for 3 
    observations to illustrate the run-wise file structure.
    
    The CTA ROOT 2D performance files contain the following histograms:
    - DiffSens_offaxis (2D) - Differential sensitivity
    - EffectiveArea_offaxis (2D) - Effective area (full containment?, reco energy?)
    - EffectiveArea80_offaxis (2D) - Effective area for 80% source containment
    - AngRes_offaxis (2D) - Angular resolution 68% containment
    - AngRes80_offaxis (2D) - Angular resolution 80% containment
    - BGRatePerSqDeg_offaxis (2D) - Background rate per square degree [counts/s/deg2]
    - BGRate_offaxis (2D) - Background rate
    - EffectiveAreaEtrue_offaxis (2D) - Effective area in true energy (finer bins)
    - MigMatrix_offaxis (3D)
    - EestOverEtrue_offaxis (3D)
    - ERes_offaxis (2D) - Energy Resolution
    - Ebias_offaxis (2D) - Energy bias
    """
    # Get database
    #entries = set_test()
    entries = set_prod1_ifae()
    #entries = set_prod2_desy()

    # Loop over entries
    for entry in entries:

        # Get attributes
        inst     = entry['inst']
        path     = entry['path']
        id       = entry['id']
        rebin    = entry['rebin']
        psftype  = entry['psftype']
        scale    = entry['scale']
        filename = path+"/"+entry['file']

        # Allocate caldb entry
        irf = caldb(inst, id)
    
        # Open calibration file
        irf.open("file")
    
        # Translate ROOT to CALDB information
        irf.root2caldb(filename, rebin=rebin, psftype=psftype, scale=scale)
