#! /usr/bin/env python
# ==========================================================================
# This script generates the TS distribution for a particular model based
# on Monte-Carlo simulations.
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
from ctools import *
from gammalib import *
import obsutils
import sys
import csv
import math


# ============== #
# cstsdist class #
# ============== #
class cstsdist(GApplication):
	"""
	This class implements the TS distribution generation script. It derives
	from the GammaLib::GApplication class which provides support for parameter
	files, command line arguments, and logging. In that way the Python
	script behaves just as a regular ctool.
	"""
	def __init__(self, *argv):
		"""
		Constructor.
		"""
		# Set name
		self.name    = "cstsdist"
		self.version = "0.1.0"
		
		# Initialise some members
		self.obs        = None
		self.bkg_model  = None
		self.full_model = None
		
		# Make sure that parfile exists
		file = self.parfile()

		# Initialise application
		if len(argv) == 0:
			GApplication.__init__(self, self.name, self.version)
		elif len(argv) ==1:
			GApplication.__init__(self, self.name, self.version, *argv)
		else:
			raise TypeError("Invalid number of arguments given.")

		# Set logger properties
		self.log_header()
		self.log.date(True)

		# Return
		return
	
	def __del__(self):
		"""
		Destructor.
		"""
		#  Write separator into logger
		if self.logTerse():
			self.log("\n")
		
		# Return
		return

	def parfile(self):
		"""
		Check if parfile exists. If parfile does not exist then create a
		default parfile. This kluge avoids shipping the cscript with a parfile.
		"""
		# Set parfile name
		parfile = self.name+".par"
		
		try:
			pars = GPars(parfile)
		except:
			# Signal if parfile was not found
			print "Parfile "+parfile+" not found. Create default parfile."
			
			# Create default parfile
			pars = GPars()
			pars.append(GPar("outfile","f","h","ts.dat","","","Output file name"))
			pars.append(GPar("ntrials","i","a","10","","","Number of trials"))
			pars.append(GPar("caldb","s","a","$GAMMALIB/share/caldb/cta","","","Calibration database"))
			pars.append(GPar("irf","s","a","kb_E_50h_v3","","","Instrument response function"))
			pars.append(GPar("type","s","a","point","","","Source model type (point/gauss/shell/disk)"))
			pars.append(GPar("index","r","h","-2.48","","","Spectral index"))
			pars.append(GPar("offset","r","a","0.0","0.0","","Source offset angle (deg)"))
			pars.append(GPar("bkg","s","a","$GAMMALIB/share/models/bkg_kb_E_50h_v3.txt","","","Background model file function (none=power law for E)"))
			pars.append(GPar("emin","r","a","0.1","0.0","","Lower energy limit (TeV)"))
			pars.append(GPar("emax","r","a","100.0","0.0","","Upper energy limit (TeV)"))
			pars.append(GPar("enumbins","i","a","0","","","Number of energy bins (0=unbinned)"))
			pars.append(GPar("duration","r","a","180000.0","","","Effective exposure time (s)"))
			pars.append(GPar("rad","r","h","5.0","","","Radius of ROI (deg)"))
			pars.append(GPar("npix","i","h","200","","","Number of pixels for binned"))
			pars.append(GPar("binsz","r","h","0.05","","","Pixel size for binned (deg/pixel)"))
			pars.append_standard()
			pars.save(parfile)
		
		# Return
		return
		
	def get_parameters(self):
		"""
		Get parameters from parfile and setup the observation.
		"""
		# Get parameters
		self.m_outfile  = self["outfile"].filename()
		self.m_ntrials  = self["ntrials"].integer()
		self.m_caldb    = self["caldb"].string()
		self.m_irf      = self["irf"].string()
		self.m_type     = self["type"].string()
		self.m_index    = self["index"].real()
		self.m_offset   = self["offset"].real()
		self.m_bkg      = self["bkg"].string()
		self.m_emin     = self["emin"].real()
		self.m_emax     = self["emax"].real()
		self.m_enumbins = self["enumbins"].integer()
		self.m_duration = self["duration"].real()
		self.m_rad      = self["rad"].real()
		self.m_npix     = self["npix"].integer()
		self.m_binsz    = self["binsz"].real()

		# Set some fixed parameters
		self.m_log   = False # Logging in client tools
		self.m_debug = False # Debugging in client tools
		
		# Setup observation
		self.obs = self.set_obs(emin=self.m_emin, emax=self.m_emax)
		
		# Initialise models. Note that we centre the point source at the Galactic
		# center as our observation is also centred at the Galactic centre, so
		# we're onaxis.
		self.bkg_model  = GModels()
		self.full_model = GModels()
		self.bkg_model.append(self.set_bkg_model())
		self.full_model.append(self.set_bkg_model())
		self.full_model.append(self.set_src_model(0.0, self.m_offset, \
		                                          flux=0.010, \
											      type=self.m_type, \
		                                          index=self.m_index))

		# Attach background model to observation container
		self.obs.models(self.bkg_model)
		
		# Return
		return
	
	def models(self, models):
		"""
		Set model.
		"""
		# Copy models
		self.model = models.copy()
	
		# Return
		return
		
	def execute(self):
		"""
		Execute the script.
		"""
		# Run the script
		self.run()
		
		# Return
		return

	def run(self):
		"""
		Run the script.
		"""
		# Switch screen logging on in debug mode
		if self.logDebug():
			self.log.cout(True)

		# Get parameters
		self.get_parameters()
		
		#  Write input parameters into logger
		if self.logTerse():
			self.log_parameters()
			self.log("\n")
		
		# Write observation into logger
		if self.logTerse():
			self.log("\n")
			self.log.header1("Observation")
			self.log(str(self.obs))
			self.log("\n")

		# Write models into logger
		if self.logTerse():
			self.log("\n")
			self.log.header1("Test model")
			self.log(str(self.full_model))
			self.log("\n")

		# Write header
		if self.logTerse():
			self.log("\n")
			self.log.header1("Generate TS distribution")

		# Loop over trials
		for seed in range(self.m_ntrials):
		
			# Make a trial
			result = self.trial(seed)
			
			# Write out result immediately
			if seed == 0:
				file = open(self.m_outfile, 'w')
				writer = csv.DictWriter(file, result['colnames'])
				writer.writerow(dict((_,_) for _ in result['colnames']))
			else:
				file = open(self.m_outfile, 'a')
			writer = csv.DictWriter(file, result['colnames'])
			writer.writerow(result['values'])
			file.close()
		
		# Return
		return
	
	def set_obs(self, lpnt=0.0, bpnt=0.0, emin=0.1, emax=100.0):
		"""
		Returns an observation container with a single CTA observation.
		
		Keywords:
		 lpnt - Galactic longitude of pointing [deg] (default: 0.0)
		 bpnt - Galactic latitude of pointing [deg] (default: 0.0)
		 emin - Minimum energy [TeV] (default: 0.1)
		 emax - Maximum energy [TeV] (default: 100.0)
		"""
		# Allocate observation container
		obs = GObservations()
	
		# Set single pointing
		pntdir = GSkyDir()
		pntdir.lb_deg(lpnt, bpnt)
		
		# Create CTA observation
		run = obsutils.set(pntdir, caldb=self.m_caldb, irf=self.m_irf, \
		                   duration=self.m_duration, \
                           emin=emin, emax=emax, rad=self.m_rad)
		
		# Append observation to container
		obs.append(run)
	
		# Return observation container
		return obs
	
	def set_bkg_model(self, fitsigma=False):
		"""
		Setup CTA background model.
		"""
		# Define radial component
		radial = GCTAModelRadialGauss(3.0)
		if fitsigma:
			radial["Sigma"].free()
		else:
			radial["Sigma"].fix()
		
		# Define spectral component
		spectrum = GModelSpectralFunc(self.m_bkg)
		
		# Create background model
		model = GCTAModelRadialAcceptance(radial, spectrum)
		model.name("Background")
		model.instruments("CTA")
	
		# Return background model
		return model
	
	def set_src_model(self, l, b, flux=1.0, index=-2.48, \
                      type="point", sigma=1.0, radius=1.0, width=0.1, \
                      fitpos=False, fitidx=False):
		"""
		Returns a single source with Crab-like spectrum. The source flux
		can be scaled in Crab units. The Crab spectrum is based on MAGIC
		observations (Albert et al. 2008, ApJ, 674, 1037).

		Parameters:
		 l      - Galactic longitude of source location [deg]
		 b      - Galactic latitude of source location [deg]
		Keywords:
		 flux   - Source flux [Crabs]
		 index  - Spectral index
		 type   - Source type ("point", "gauss", "disk", "shell")
		 sigma  - Gaussian sigma (for type="gauss")
		 radius - Disk or shell inner radius [deg] (for type="disk" and type="shell")
		 width  - Shell width [deg] (for type="shell")
		 fitpos - Fit position? (default: True)
		 fitidx - Fit index? (default: True)
		"""
		# Set source location
		location = GSkyDir()
		location.lb_deg(l, b)
	
		# Set source spectrum
		spectrum = GModelSpectralPlaw(flux, index)
		spectrum["Prefactor"].scale(5.7e-16)
		spectrum["PivotEnergy"].value(0.3)
		spectrum["PivotEnergy"].scale(1.0e6)
		if fitidx:
			spectrum["Index"].free()
		else:
			spectrum["Index"].fix()	

		# Set source
		if type == "point":
			spatial = GModelSpatialPtsrc(location)
			if fitpos:
				spatial[0].free()
				spatial[1].free()
			source  = GModelPointSource(spatial, spectrum)
		elif type == "gauss":
			radial = GModelRadialGauss(location, sigma)
			if fitpos:
				radial[0].free()
				radial[1].free()
			source = GModelExtendedSource(radial, spectrum)
		elif type == "disk":
			radial = GModelRadialDisk(location, radius)
			if fitpos:
				radial[0].free()
				radial[1].free()
			source = GModelExtendedSource(radial, spectrum)
		elif type == "shell":
			radial = GModelRadialShell(location, radius, width)
			if fitpos:
				radial[0].free()
				radial[1].free()
			source = GModelExtendedSource(radial, spectrum)
		else:
			self.log("ERROR: Unknown source type '"+type+"'.\n")
			return None
	
		# Set source name
		source.name("Test")
	
		# Return source
		return source
	
	def trial(self, seed):
		"""
		Create the TS for a single trial.
		
		Parameters:
		 seed - Random number generator seed
		"""
		# Write header
		if self.logExplicit():
			self.log.header2("Trial "+str(seed+1))

		# Simulate events
		sim = obsutils.sim(self.obs, \
		                   nbins=self.m_enumbins, \
		                   seed=seed, \
		                   binsz=self.m_binsz, \
		                   npix=self.m_npix, \
		                   log=self.m_log, debug=self.m_debug)

		# Determine number of events in simulation
		nevents = 0.0
		for run in sim:
			nevents += run.events().number()

		# Write simulation results
		if self.logExplicit():
			self.log.header3("Simulation")
			self.log.parformat("Number of simulated events")
			self.log(nevents)
			self.log("\n")

		# Fit background only
		sim.models(self.bkg_model)
		like_bgm   = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
		result_bgm = like_bgm.obs().models()
		LogL_bgm   = like_bgm.opt().value()
		npred_bgm  = like_bgm.obs().npred()

		# Write background fit results
		if self.logExplicit():
			self.log.header3("Background model fit")
			self.log.parformat("log likelihood")
			self.log(LogL_bgm)
			self.log("\n")
			self.log.parformat("Number of predicted events")
			self.log(npred_bgm)
			self.log("\n")
			for model in result_bgm:
				self.log.parformat("Model")
				self.log(model.name())
				self.log("\n")
				for par in model:
					self.log(str(par)+"\n")

		# Fit background and test source
		sim.models(self.full_model)
		like_all   = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
		result_all = like_all.obs().models()
		LogL_all   = like_all.opt().value()
		npred_all  = like_all.obs().npred()
		ts         = 2.0*(LogL_bgm-LogL_all)
		
		# Write background and test source fit results
		if self.logExplicit():
			self.log.header3("Background and test source model fit")
			self.log.parformat("Test statistics")
			self.log(ts)
			self.log("\n")
			self.log.parformat("log likelihood")
			self.log(LogL_all)
			self.log("\n")
			self.log.parformat("Number of predicted events")
			self.log(npred_all)
			self.log("\n")
			for model in result_all:
				self.log.parformat("Model")
				self.log(model.name())
				self.log("\n")
				for par in model:
					self.log(str(par)+"\n")
					
		# Write result
		elif self.logTerse():
			self.log.parformat("Trial "+str(seed))
			self.log("TS=")
			self.log(ts)
			self.log("  Prefactor=")
			self.log(result_all["Test"]["Prefactor"].value())
			self.log("+/-")
			self.log(result_all["Test"]["Prefactor"].error())
			self.log("\n")
		
		# Initialise results
		colnames = []
		values   = {}
		
		# Set TS value
		colnames.append("TS")
		values["TS"] = ts

		# Set logL for background fit
		colnames.append("LogL_bgm")
		values["LogL_bgm"] = LogL_bgm

		# Set logL for full fit
		colnames.append("LogL_all")
		values["LogL_all"] = LogL_all

		# Set Nevents
		colnames.append("Nevents")
		values["Nevents"] = nevents

		# Set Npred for background fit
		colnames.append("Npred_bkg")
		values["Npred_bkg"] = npred_bgm

		# Set Npred for full fit
		colnames.append("Npred_all")
		values["Npred_all"] = npred_all
		
		# Gather free full fit parameters
		for i in range(result_all.size()):
			model      = result_all[i]
			model_name = model.name()
			for k in range(model.size()):
				par = model[k]
				if par.isfree():
				
					# Set parameter name
					name = model_name+"_"+par.name()
					
					# Append value
					colnames.append(name)
					values[name] = par.value()
					
					# Append error
					name = "Unc_"+name
					colnames.append(name)
					values[name] = par.error()
		
		# Bundle together results
		result = {'colnames': colnames, 'values': values}
		
		# Return
		return result


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
	"""
	Generates TS distribution.
	"""
	# Create instance of application
	app = cstsdist(sys.argv)
	
	# Open logfile
	app.logFileOpen()
	
	# Execute application
	app.execute()
	