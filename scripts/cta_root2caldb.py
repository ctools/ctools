#! /usr/bin/env python
# ==========================================================================
# This script generates IRFs in the CALDB format of HEASARC using ROOT 2D
# performance files.
#
# Requirements:
# - pyroot
#
# Copyright (C) 2011 Jurgen Knodlseder
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
	def __init__(self, conf, path=""):
		"""
		Constructor.
		"""
		# Store root path
		self.path = path
		
		# Initialise some members
		self.cif     = None
		self.ea      = None
		self.psf     = None
		self.edisp   = None
		self.hdu_cif = None

		# Store UTC date of job execution
		self.utc = datetime.utcnow().isoformat()[:19]
		
		# Set CALDB information
		self.cal_tel      = "CTA"
		self.cal_inst     = conf.upper()
		self.cal_det      = "NONE"
		self.cal_flt      = "NONE"
		self.cal_class    = "BCF"
		self.cal_type     = "DATA"
		self.cal_qual     = 0            # 0=good, 1=bad, 2=dubious, ...
		self.cal_date     = "11/01/01"
		self.val_date     = "2011-01-01"
		self.val_time     = "00:00:00"
		self.ref_time     = 51544.0
		self.cal_version  = "VERSION(DC1)"
		self.cal_cut      = "CLASS(BEST)"
		self.cal_analysis = "ANALYSIS(IFAE)"
		#
		self.ea_name      = "EFF_AREA"
		self.ea_doc       = "CAL/GEN/92-019"
		self.ea_bounds    = [self.cal_version, self.cal_cut, self.cal_analysis]
		self.ea_desc      = "CTA effective area"
		#
		self.psf_name     = "RPSF"
		self.psf_doc      = "CAL/GEN/92-020"
		self.psf_bounds   = [self.cal_version, self.cal_cut, self.cal_analysis]
		self.psf_desc     = "CTA point spread function"
		#
		self.edisp_name   = "EDISP"
		self.edisp_doc    = "???"
		self.edisp_bounds = [self.cal_version, self.cal_cut, self.cal_analysis]
		self.edisp_desc   = "CTA energy dispersion"
			
		# Create directory structure
		self.make_dirs()
		
		# Return
		return
	
	def __del__(self):
		"""
		"""
		# Close calibration
		self.close()

		# Return
		return
	
	def make_dirs(self):
		"""
		"""
		# Set base path
		self.base_dir  = "data"
		self.base_dir += "/"+self.cal_tel.lower()
		self.base_dir += "/"+self.cal_inst.lower()


		# Set directory names
		self.bcf_dir   = self.base_dir+"/bcf"
		self.ea_dir    = self.bcf_dir+"/ea"
		self.psf_dir   = self.bcf_dir+"/psf"
		self.edisp_dir = self.bcf_dir+"/edisp"
		
		# Add path (if required)
		if len(self.path) > 0:
			self.base_path  = self.path+"/"+self.base_dir
			self.bcf_path   = self.path+"/"+self.bcf_dir
			self.ea_path    = self.path+"/"+self.ea_dir
			self.psf_path   = self.path+"/"+self.psf_dir
			self.edisp_path = self.path+"/"+self.edisp_dir
		else:
			self.base_path  = self.base_dir
			self.bcf_path   = self.bcf_dir
			self.ea_path    = self.ea_dir
			self.psf_path   = self.psf_dir
			self.edisp_path = self.edisp_dir

		# Create directory structure
		if not os.path.isdir(self.ea_path):
			os.makedirs(self.ea_path)
		if not os.path.isdir(self.psf_path):
			os.makedirs(self.psf_path)
		if not os.path.isdir(self.edisp_path):
			os.makedirs(self.edisp_path)
	
		# Return
		return
	
	def open(self, name, clobber=True):
		"""
		Open existing or create new calibration.
		"""
		# Set calibrate file names
		self.ea_file    = "ea_"+name+".fits"
		self.psf_file   = "psf_"+name+".fits"
		self.edisp_file = "edisp_"+name+".fits"
		
		# Open calibration database index
		self.cif = GFits(self.base_path+"/caldb.indx", True)
		
		# If file has no CIF extension than create it now
		try:
			self.hdu_cif = self.cif.table("CIF")
		except:
			self.create_cif()
			self.hdu_cif = self.cif.table("CIF")
		
		# Open calibration files
		self.open_ea(self.ea_path+"/"+self.ea_file)
		self.open_psf(self.psf_path+"/"+self.psf_file)
		self.open_edisp(self.edisp_path+"/"+self.edisp_file)
		
		# Return
		return

	def close(self):
		"""
		Close calibration.
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
		
		# Close effective area
		if self.ea != None:
			self.ea.save(True)
			self.ea.close()
			self.ea = None

		# Close point spread function
		if self.psf != None:
			self.psf.save(True)
			self.psf.close()
			self.psf = None

		# Close energy dispersion
		if self.edisp != None:
			self.edisp.save(True)
			self.edisp.close()
			self.edisp = None
		
		# Return
		return
	
	def create_cif(self):
		"""
		Create CIF extension in FITS file.
		"""
		# Create binary table
		table = GFitsBinTable()
		
		# Set boundary
		
		# Attach columns. Reference: CAL/GEN/92-008
		table.append_column(GFitsTableStringCol("TELESCOP", 0, 10))
		table.append_column(GFitsTableStringCol("INSTRUME", 0, 10))
		table.append_column(GFitsTableStringCol("DETNAM", 0, 20))
		table.append_column(GFitsTableStringCol("FILTER", 0, 10))
		table.append_column(GFitsTableStringCol("CAL_DEV", 0, 20))
		table.append_column(GFitsTableStringCol("CAL_DIR", 0, 70))
		table.append_column(GFitsTableStringCol("CAL_FILE", 0, 40))
		table.append_column(GFitsTableStringCol("CAL_CLAS", 0, 3))
		table.append_column(GFitsTableStringCol("CAL_DTYP", 0, 4))
		table.append_column(GFitsTableStringCol("CAL_CNAM", 0, 20))
		table.append_column(GFitsTableStringCol("CAL_CBD", 0, 70, 9))
		table.append_column(GFitsTableShortCol("CAL_XNO", 0))
		table.append_column(GFitsTableStringCol("CAL_VSD", 0, 10))
		table.append_column(GFitsTableStringCol("CAL_VST", 0, 8))
		table.append_column(GFitsTableDoubleCol("REF_TIME", 0))
		table.append_column(GFitsTableShortCol("CAL_QUAL", 0))
		table.append_column(GFitsTableStringCol("CAL_DATE", 0, 8))
		table.append_column(GFitsTableStringCol("CAL_DESC", 0, 70))
		
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
		"""
		# Append 3 rows to CIF extension
		self.hdu_cif.append_rows(3)
		
		# Add generic information for these 3 rows
		for i in range(3):
		
			# Set row number
			row = i+self.hdu_cif.nrows()-3
			
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
		row = self.hdu_cif.nrows()-3
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
		row = self.hdu_cif.nrows()-2
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
		row = self.hdu_cif.nrows()-1
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
		
		# Return
		return

	def open_ea(self, filename):
		"""
		Open effective area extension. If no extension was found
		then create one. We do not yet set any data as we don't
		know the dimensions of the table yet.
		"""
		# Open FITS file
		self.ea = GFits(filename, True)

		# Get extensions. If they do not exist then create them now
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

	def open_psf(self, filename):
		"""
		Open point spread function extension. If no extension was
		found then create one. We do not yet set any data as we don't
		know the dimensions of the table yet.
		"""
		# Open FITS file
		self.psf = GFits(filename, True)

		# Get extensions. If they do not exist then create them now
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

	def open_edisp(self, filename):
		"""
		Open energy dispersion extension. If no extension was found
		then create one. We do not yet set any data as we don't
		know the dimensions of the table yet.
		"""
		# Open FITS file
		self.edisp = GFits(filename, True)

		# Get extensions. If they do not exist then create them now
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

	def set_ogip_keywords(self, hdu, hdudoc, hduclas):
		"""
		Set standard OGIP keywords for extension.
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
		ROOT 2D histogram.
		
		Parameters:
		 array - ROOT 2D histgram.
		 hdu   - FITS HDU.
		Keyword:
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
		if not hdu.hascolumn("ENERG_LO"):
			hdu.append_column(GFitsTableFloatCol("ENERG_LO", 1, neng))
			hdu["ENERG_LO"].unit("TeV")
			for ieng in range(neng):
				e_lo = pow(10.0, energies.GetBinLowEdge(ieng+1))
				hdu["ENERG_LO"][0,ieng] = e_lo
		#
		# ENERG_HI
		if not hdu.hascolumn("ENERG_HI"):
			hdu.append_column(GFitsTableFloatCol("ENERG_HI", 1, neng))
			hdu["ENERG_HI"].unit("TeV")
			for ieng in range(neng):
				e_hi = pow(10.0, energies.GetBinUpEdge(ieng+1))
				hdu["ENERG_HI"][0,ieng] = e_hi
		#
		# THETA_LO
		if not hdu.hascolumn("THETA_LO"):
			hdu.append_column(GFitsTableFloatCol("THETA_LO", 1, noffset))
			hdu["THETA_LO"].unit("deg")
			for ioff in range(noffset):
				o_lo = offsets.GetBinLowEdge(ioff+1)
				hdu["THETA_LO"][0,ioff] = o_lo

		#
		# THETA_LO
		if not hdu.hascolumn("THETA_HI"):
			hdu.append_column(GFitsTableFloatCol("THETA_HI", 1, noffset))
			hdu["THETA_HI"].unit("deg")
			for ioff in range(noffset):
				o_hi = offsets.GetBinUpEdge(ioff+1)
				hdu["THETA_HI"][0,ioff] = o_hi
		#
		# "NAME"
		if not hdu.hascolumn(name):
			hdu.append_column(GFitsTableFloatCol(name, 1, neng*noffset))
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
	
	def root2caldb(self, filename):
		"""
		Translate ROOT to CALDB information.
		
		Parameters:
		 filename - ROOT 2D performance filename.
		"""
		# Open performance file
		file = TFile(filename)

		# Allocate ROOT 2D array
		array = TH2F()

		# Make 2D effective 80% area
		if self.hdu_ea != None:
			file.GetObject("EffectiveArea80_offaxis", array)
			bounds = self.make_2D(array, self.hdu_ea, "EFFAREA", "m2", scale=1.0/0.8)
			for b in bounds:
				self.ea_bounds.append(b)
			self.set_cif_keywords(self.hdu_ea, self.ea_name, \
			                      self.ea_bounds, self.ea_desc)

		# Make 2D 68% point spread function
		if self.hdu_psf != None:
			file.GetObject("AngRes_offaxis", array)
			bounds = self.make_2D(array, self.hdu_psf, "R68", "deg")
			for b in bounds:
				self.psf_bounds.append(b)
			self.set_cif_keywords(self.hdu_psf, self.psf_name, \
			                      self.psf_bounds, self.psf_desc)

		# Make 2D 80% point spread function
		if self.hdu_psf != None:
			file.GetObject("AngRes80_offaxis", array)
			self.make_2D(array, self.hdu_psf, "R80", "deg")

		# Make 2D energy resolution
		if self.hdu_edisp != None:
			file.GetObject("ERes", array)
			bounds = self.make_2D(array, self.hdu_edisp, "ERES68", "")
			for b in bounds:
				self.edisp_bounds.append(b)
			self.set_cif_keywords(self.hdu_edisp, self.edisp_name, \
			                      self.edisp_bounds, self.edisp_desc)
		
		# Return
		return

	def set_cif_keywords(self, hdu, name, bounds,desc):
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
		

#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	This script translates a ROOT 2D performance file in a CALDB compliant
	CTA response.
	"""
	# Set ROOT filename
	#path     = "/Users/jurgen/Documents/Travail/Projects/CTA/WP-MC/root/IFAEOffaxis"
	#irf      = "SubarrayE_offaxis.root"
	path     = "/Users/jurgen/Documents/Travail/Projects/CTA/WP-MC/root/IFAEOffaxisPerformanceBEI_Nov2011"
	irf      = "SubarrayE_IFAE_50hours_20111121_offaxis.root"
	filename = path+"/"+irf

	# Allocate caldb
	db = caldb("E")
	
	# Open calibration
	db.open("test")
	
	# Translate ROOT to CALDB information
	db.root2caldb(filename)
	