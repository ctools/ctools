#! /usr/bin/env python
# ==========================================================================
# This script generates CALDB compliant performance tables from MC ROOT
# files.
#
# Note that the effective area given in the ASCII file is the full effective
# area (i.e. the area before applying any theta cut). This is required by
# GammaLib as the maximum likelihood analysis does not imply any theta cut.
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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
#from ROOT import TFile, TH1F, TH2F
from ROOT import TFile, TH1F
from datetime import datetime
import gammalib
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
    def __init__(self, inst, path=""):
        """
        CALDB constructor. The constructor will create one CALDB entry for
        an observation.

        Parameters:
         inst  - Instrument string.
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
        self.cal_det      = "NONE"
        self.cal_flt      = "NONE"
        self.cal_class    = "BCF"
        self.cal_type     = "DATA"
        self.cal_qual     = 0            # 0=good, 1=bad, 2=dubious, ...
        self.cal_date     = "14/01/01"
        self.val_date     = "2014-01-01"
        self.val_time     = "00:00:00"
        self.ref_time     = 51544.0
        self.cal_name     = "NAME()"
        self.cal_version  = "VERSION(Prod1)"
        self.cal_cut      = "CLASS(BEST)"
        self.cal_analysis = "ANALYSIS(MPIK)"
        #
        self.ea_name      = "EFF_AREA"
        self.ea_doc       = "???"
        self.ea_bounds    = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]
        self.ea_desc      = "CTA effective area"
        #
        self.psf_name     = "RPSF"
        self.psf_doc      = "???"
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

            data/<tel>/<inst>/bcf/

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

        # Optionally prefix path
        if len(self.path) > 0:
            self.base_path  = self.path+"/"+self.base_dir
            self.bcf_path   = self.path+"/"+self.bcf_dir
        else:
            self.base_path  = self.base_dir
            self.bcf_path   = self.bcf_dir

        # Create directory structure
        if not os.path.isdir(self.bcf_path):
            os.makedirs(self.bcf_path)

        # Return
        return

    def add(self, rspname, split=False, clobber=True):
        """
        Add new calibration. The actual version will put
        all calibrations in the same file, although each part of the response
        function will have its own logical name. We can thus easily modify
        the script to put each calibration information in a separate file.

        Parameters:
         rspname - Response name
        Keywords:
         split   - Split IRF over several files?
         clobber - Overwrite existing files?
        """
        # Set calibrate file names
        if split:
            self.ea_file    = "ea_"+rspname+".dat"
            self.psf_file   = "psf_"+rspname+".dat"
            self.edisp_file = "edisp_"+rspname+".dat"
            self.bgd_file   = "bgd_"+rspname+".dat"
        else:
            self.ea_file    = "irf_"+rspname+".dat"
            self.psf_file   = "irf_"+rspname+".dat"
            self.edisp_file = "irf_"+rspname+".dat"
            self.bgd_file   = "irf_"+rspname+".dat"

        # Open calibration database index
        if self.cif == None:
            self.cif = gammalib.GFits(self.base_path+"/caldb.indx", True)

        # If file has no CIF extension than create it now
        try:
            self.hdu_cif = self.cif.table("CIF")
        except:
            self.create_cif()
            self.hdu_cif = self.cif.table("CIF")

        # Set response name
        self.cal_name = "NAME("+rspname+")"

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
        #self.add_cif_info()

        # Close CIF
        if self.cif != None:
            self.cif.save(True)
            self.cif.close()
            self.cif = None

        # Free
        self.cif = None

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
        # Set boundaries
        self.ea_bounds    = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]
        self.psf_bounds   = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]
        self.edisp_bounds = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]
        self.bgd_bounds   = [self.cal_name, self.cal_version, self.cal_cut, self.cal_analysis]

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
        self.hdu_cif["CAL_DIR"][row]   = self.bcf_dir
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
        self.hdu_cif["CAL_DIR"][row]   = self.bcf_dir
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
        self.hdu_cif["CAL_DIR"][row]   = self.bcf_dir
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
        self.hdu_cif["CAL_DIR"][row]   = self.bcf_dir
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

    def root2caldb(self, filename):
        """
        Generate files.

        Parameters:
         filename - Performance file name
        """
        # Open ROOT performance file
        rootfile = TFile(filename)

        # Set filenames
        ea_filename    = self.bcf_path+"/"+self.ea_file
        psf_filename   = self.bcf_path+"/"+self.psf_file
        edisp_filename = self.bcf_path+"/"+self.edisp_file
        bgd_filename   = self.bcf_path+"/"+self.bgd_file

        # Open ASCII performance tables
        if self.ea_file == self.psf_file and \
           self.ea_file == self.edisp_file and \
           self.ea_file == self.bgd_file :
            split = False
            fprf  = open(ea_filename, "w")
            files = [fprf]
        else:
            split = True
            fea    = open(ea_filename, "w")
            fpsf   = open(psf_filename, "w")
            fedisp = open(edisp_filename, "w")
            fbgd   = open(bgd_filename, "w")
            files  = [fea, fpsf, fedisp, fbgd]

        # Write header
        for file in files:
            file.write("log(E) Area  r68  r80  ERes. BG Rate  Diff Sens\n")

        # Get histograms. Konrad's root files do not have an EffectiveArea80
        # histogram, but his effective area is for 80% containment radius
        # anyways, so we simply read this histogram into aeff80.
        sens        = TH1F()
        bgrate      = TH1F()
        bgratesqdeg = TH1F()
        aeff        = TH1F()
        aeff80      = TH1F()
        angres      = TH1F()
        angres80    = TH1F()
        eres        = TH1F()
        rootfile.GetObject("DiffSens", sens)
        rootfile.GetObject("BGRate", bgrate)
        rootfile.GetObject("BGRatePerSqDeg", bgratesqdeg)
        rootfile.GetObject("EffectiveArea", aeff)
        try:
            rootfile.GetObject("EffectiveArea80", aeff80)
        except:
            rootfile.GetObject("EffectiveArea", aeff80)
        rootfile.GetObject("AngRes", angres)
        rootfile.GetObject("AngRes80", angres80)
        rootfile.GetObject("ERes", eres)

        # Setup vectors
        energies  = bgrate.GetXaxis()
        nbins_eng = energies.GetNbins()
        for i in range(nbins_eng):

            # Get logE and background rate per squared degree
            logE   = energies.GetBinCenter(i+1)
            bgd    = bgrate.GetBinContent(i+1,1)
            bgdsq  = bgratesqdeg.GetBinContent(i+1,1)
            area   = aeff.GetBinContent(i+1,1)
            area80 = aeff80.GetBinContent(i+1,1)
            r68    = angres.GetBinContent(i+1,1)
            r80    = angres80.GetBinContent(i+1,1)
            dE     = eres.GetBinContent(i+1,1)
            diffs  = sens.GetBinContent(i+1,1)

            # Skip any NaNs
            if math.isnan(area80):
                continue

            # Compute energy (in MeV)
            energy = math.pow(10.0, logE)*1.0e6
            emin   = math.pow(10.0, logE-0.1)*1.0e6
            emax   = math.pow(10.0, logE+0.1)*1.0e6
            ewidth = emax-emin

            # Convert into background rate per steradian and MeV
            bkg_rate  = bgdsq / (0.01745329*0.01745329) / ewidth

            # Compute control background rate (this only works for KB files
            # as they are for 80% containment radius). But we don't use this
            # anyways, it's just for control on display ...
            omega = 2.0 * math.pi * (1.0 - math.cos(math.radians(r80)))
            if omega > 0.0:
                bkg_rate2 = bgd / omega / ewidth
            else:
                bkg_rate2 = 0.0

            # Compute full effective area by dividing area80/0.80
            area = area80/0.80

            # Write results in file
            line = "%.1f  %.1f  %.4f  %.4f %.4f %.5e  %.5e" % \
                   (logE, area, r68, r80, dE, bgd, diffs)
            for file in files:
                file.write(line+"\n")
            #fbkg.write(str(energy)+" "+str(bkg_rate)+"\n")

        # Write trailer
        for file in files:
            file.write("---------------------------------------------\n")
            file.write("Notes\n")
            file.write(" 1) log(E) = log10(E/TeV) - bin centre\n")
            file.write(" 2) Eff Area - in square metres after background cut (no theta cut)\n")
            file.write(" 3) Ang. Res - 68% containment radius of gamma-ray PSF post cuts - in degrees\n")
            file.write(" 4) Ang. Res - 80% containment radius of gamma-ray PSF post cuts - in degrees\n")
            file.write(" 5) Fractional Energy Resolution (rms)\n")
            file.write(" 6) BG Rate  - inside point-source selection region - post call cuts - in Hz\n")
            file.write(" 7) Diff Sens - differential sensitivity for this bin expressed as E^2 dN/dE\n")
            file.write(" - in erg cm^-2 s^-1 - for a 50 hours exposure - 5 sigma significance including\n")
            file.write(" systematics and statistics and at least 10 photons.\n")

        # Close files
        for file in files:
            file.close()

        # Close ROOT performance file
        rootfile.Close()

        # Add information to CIF. We do this now as before we did not
        # necessarily have all the information at hand (in particular
        # about the boundaries)
        self.add_cif_info()

        # Return
        return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    """
    # Get ROOT filename from command line argument
    #if (len(sys.argv) > 1):
    #    filename = sys.argv[1]
    #else:
    db = "/Users/jurgen/Documents/Travail/Projects/CTA/WP-MC/root/Kb_all_20deg_root"
    #db = "/Users/jurgen/Documents/Travail/Projects/CTA/WP-MC/root/Kb_all_50deg_root"

    # Get list of all performance files
    files = glob.glob(db+"/*.root")
    files.sort()

    # Allocate caldb
    irf = caldb("kb")

    # Loop over all files
    for f in files:

        # Get filename
        head, tail = os.path.split(f)

        # Build name
        name = tail.strip(".root")
        if name.find("kb_") != -1:
            name = name.replace("_20deg", "")
        elif name.find("IFAE_") != -1:
            name  = name.replace("Subarray", "")
            name  = name.replace("IFAE_", "")
            name  = "IFAE_"+name
            index = name.find("hours")
            if index != -1:
                name = name[0:index+1]
        print(tail)

        # Open calibration file
        irf.add(name)

        # Translate ROOT to CALDB information
        irf.root2caldb(f)
