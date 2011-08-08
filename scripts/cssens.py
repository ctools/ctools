#! /usr/bin/env python
# ==========================================================================
# This script computes the CTA sensitivity for a number of energy bins using
# ctlike. Crab spectra are fitted in narrow energy bins to simulated data,
# and the flux level is determined that leads to a particular significance.
#
# The significance can be either determined using the Test Statistic (which
# is defined as twice the likelihood difference between fitting with and
# without the test source) or using the source model scaling factor (or
# Prefactor) divided by its statistical error. The latter has been turned
# out to be more reliable, and is used as the default parameter.
#
# As background model, a GCTAModelRadialGauss model is used for the spatial
# component. For the spectral component, either a file function is used
# (specified by the "bkg" parameter), or if left blank, a power law adjusted
# to Konrad's configuration "E" performance file is employed.
#
# The source is modelled as a point source
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
from ctatools import *
from gammalib import *
import obsutils
import sys
import csv
import math


# ============ #
# cssens class #
# ============ #
class cssens(GApplication):
	"""
	This class implements the sensitivity computation script. It derives
	from the GammaLib::GApplication class which provides support for parameter
	files, command line arguments, and logging. In that way the Python
	script behaves just as a regular ctool.
	"""
	def __init__(self, *argv):
		"""
		Constructor.
		"""
		# Set name
		self.name    = "cssens"
		self.version = "0.3.0"
		
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
		self.logFileOpen()
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
			pars.append(GPar("outfile","f","h","sensitivity.dat","","","Output file name"))
			pars.append(GPar("caldb","s","a","$GAMMALIB/share/caldb/cta","","","Calibration database"))
			pars.append(GPar("irf","s","a","kb_E_50h_v3","","","Instrument response function"))
			pars.append(GPar("type","s","a","point","","","Source model type (point/gauss/shell/disk)"))
			pars.append(GPar("offset","r","a","0.0","0.0","","Source offset angle (deg)"))
			pars.append(GPar("bkg","s","a","$GAMMALIB/share/models/bkg_kb_E_50h_v3.txt","","","Background model file function (none=power law for E)"))
			pars.append(GPar("duration","r","a","180000.0","","","Effective exposure time (s)"))
			pars.append(GPar("rad","r","a","5.0","","","Radius of ROI (deg)"))
			pars.append(GPar("enumbins","i","h","0","","","Number of energy bins (0=unbinned)"))
			pars.append(GPar("npix","i","h","200","","","Number of pixels for binned"))
			pars.append(GPar("binsz","r","h","0.05","","","Pixel size for binned (deg/pixel)"))
			pars.append(GPar("sigma","r","h","5.0","","","Significance threshold"))
			pars.append(GPar("ts_use","b","h","yes","","","Use TS?"))
			pars.append(GPar("index","r","h","-2.48","","","Assumed spectral index"))
			pars.append(GPar("radius","r","h","0.1","","","Extended source model radius"))
			pars.append(GPar("width","r","h","0.05","","","Extended source model width"))
			pars.append(GPar("num_avg","i","h","3","","","Number of iterations for sliding average"))
			pars.append_standard()
			pars.save(parfile)
		
		# Return
		return
		
	def get_parameters(self):
		"""
		Get parameters from parfile.
		"""
		# Get parameters
		self.m_outfile  = self["outfile"].filename()
		self.m_caldb    = self["caldb"].string()
		self.m_irf      = self["irf"].string()
		self.m_type     = self["type"].string()
		self.m_offset   = self["offset"].real()
		self.m_bkg      = self["bkg"].string()
		self.m_duration = self["duration"].real()
		self.m_roi      = self["rad"].real()
		self.m_nbins    = self["enumbins"].integer()
		self.m_npix     = self["npix"].integer()
		self.m_binsz    = self["binsz"].real()
		self.m_ts_thres = self["sigma"].real()*self["sigma"].real()
		self.m_use_ts   = self["ts_use"].boolean()
		self.m_index    = self["index"].real()
		self.m_radius   = self["radius"].real()
		self.m_width    = self["width"].real()
		self.m_max_iter = self["max_iter"].integer()
		self.m_num_avg  = self["num_avg"].integer()

		# Set some fixed parameters
		self.m_log   = False # Logging in client tools
		self.m_debug = False # Debugging in client tools
		
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
		
		# Initialise script
		colnames = ['loge', 'emin', 'emax', 'crab', 'pflux', 'eflux', 'diffSens']
		results  = []

		# Initialise models. Note that we centre the point source at the Galactic
		# center as our observation is also centred at the Galactic centre, so
		# we're onaxis.
		bkg_model  = GModels()
		full_model = GModels()
		bkg_model.append(self.set_bkg_model())
		full_model.append(self.set_bkg_model())
		full_model.append(self.set_src_model(0.0, \
		                                     self.m_offset, \
											 type=self.m_type, \
		                                     index=self.m_index))

		# Write models into logger
		if self.logTerse():
			self.log("\n")
			self.log.header1("Models")
			self.log.header2("Background model")
			self.log(str(bkg_model))
			self.log("\n\n")
			self.log.header2("Full model")
			self.log(str(full_model))
			self.log("\n")

		# Write heder
		if self.logTerse():
			self.log("\n")
			self.log.header1("Sensitivity determination")

		# Loop over energy bands. The number of energy band is still
		# fixed ...
		for ieng in range(21):
		
			# Set energies
			loge  = -1.7 + ieng * 0.2
			emean = pow(10.0, loge)
			emin  = pow(10.0, loge-0.1)
			emax  = pow(10.0, loge+0.1)

			# Setup observation(s) 
			obs = self.set_obs(emin=emin, emax=emax)
			
			# Determine sensitivity
			result = self.get_sensitivity(obs, bkg_model, full_model)
			
			# Write results
			if ieng == 0:
				f      = open(self.m_outfile, 'w')
				writer = csv.DictWriter(f, colnames)
				writer.writerow(dict((_,_) for _ in colnames))
				writer.writerow(result)
				f.close()
			else:
				f = open(self.m_outfile, 'a')
				writer = csv.DictWriter(f, colnames)
				writer.writerow(result)
				f.close()
			
			# Append results
			results.append(result)
		
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
                           emin=emin, emax=emax, rad=self.m_roi)
		
		# Append observation to container
		obs.append(run)
	
		# Return observation container
		return obs
	
	def set_bkg_model(self, fitidx=False):
		"""
		Setup standard CTA background model. We use a simple power law here,
		scaled to Konrad's E configuration performance table. 
		
		Keywords:
		 fitidx - Fit spectral index (default: False)
		"""
		# Define radial component
		radial = GCTAModelRadialGauss(3.0)
		
		# Define spectral component. If a background model is given then
		# use a file function. Otherwise use a static power law.
		if len(self.m_bkg) > 0:
			spectrum = GModelSpectralFunc(self.m_bkg)
		else:
			spectrum = GModelSpectralPlaw(61.8, -1.85)
			spectrum["Prefactor"].scale(1.0e-6)
			spectrum["PivotEnergy"].value(1.0)
			spectrum["PivotEnergy"].scale(1.0e6)
			if fitidx:
				spectrum["Index"].free()
			else:
				spectrum["Index"].fix()
		
		# Create background model
		model = GCTAModelRadialAcceptance(radial, spectrum)
		model.name("Background")
		model.instruments("CTA")
	
		# Return background model
		return model
	
	def set_src_model(self, l, b, flux=1.0, index=-2.48, \
                      type="point", fitpos=False, fitidx=False):
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
		 fitpos - Fit position and size? (default: False)
		 fitidx - Fit index? (default: False)
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
				spatial["RA"].free()
				spatial["DEC"].free()
			else:
				spatial["RA"].fix()
				spatial["DEC"].fix()
			source  = GModelPointSource(spatial, spectrum)
		elif type == "gauss":
			radial = GModelRadialGauss(location, self.radius)
			if fitpos:
				radial["RA"].free()
				radial["DEC"].free()
				radial["Sigma"].free()
			else:
				radial["RA"].fix()
				radial["DEC"].fix()
				radial["Sigma"].fix()
			source = GModelExtendedSource(radial, spectrum)
		elif type == "disk":
			radial = GModelRadialDisk(location, self.radius)
			if fitpos:
				radial["RA"].free()
				radial["DEC"].free()
				radial["Radius"].free()
			else:
				radial["RA"].fix()
				radial["DEC"].fix()
				radial["Radius"].fix()
			source = GModelExtendedSource(radial, spectrum)
		elif type == "shell":
			radial = GModelRadialShell(location, self.radius,  self.width)
			if fitpos:
				radial["RA"].free()
				radial["DEC"].free()
				radial["Radius"].free()
				radial["Width"].free()
			else:
				radial["RA"].fix()
				radial["DEC"].fix()
				radial["Radius"].fix()
				radial["Width"].fix()
			source = GModelExtendedSource(radial, spectrum)
		else:
			self.log("ERROR: Unknown source type '"+type+"'.\n")
			return None
	
		# Set source name
		source.name("Test")
	
		# Return source
		return source
	
	def get_sensitivity(self, obs, bkg_model, full_model):
		"""
		Determine sensitivity for a given observations.
		
		Parameters:
		 obs        - Observation container
		 bkg_model  - Background model
		 full_model - Source model
		"""
		# Set TeV->erg conversion factor
		tev2erg = 1.6021764
		
		# Determine energy boundaries from first observation in
		# the container
		for run in obs:
			emin      = run.events().ebounds().emin()
			emax      = run.events().ebounds().emax()
			loge      = math.log10(math.sqrt(emin.TeV()*emax.TeV()))
			erg_mean  = math.pow(10.0, loge) * tev2erg
			erg_width = (emax.TeV()-emin.TeV()) * tev2erg
			break

		# Write header
		if self.logTerse():
			self.log("\n")
			self.log.header2("Energies: "+str(emin)+" - "+str(emax))

		# Initialise loop
		flux_value     = []
		pflux_value    = []
		eflux_value    = []
		diffSens_value = []
		iter           = 0
		test_flux      = 0.1  # This is the initial test flux in Crab units
		
		# Loop until we break
		while True:
		
			# Update iteration counter
			iter += 1
			
			# Write header
			if self.logExplicit():
				self.log.header2("Iteration "+str(iter))

			# Set test flux
			full_model["Test"]['Prefactor'].value(test_flux)
			obs.models(full_model)
			
			# Simulate events
			sim = obsutils.sim(obs, nbins=self.m_nbins, seed=iter, \
			                   binsz=self.m_binsz, npix=self.m_npix, \
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
			sim.models(bkg_model)
			like       = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
			result_bgm = like.obs().models().copy()
			LogL_bgm   = like.opt().value()
			npred_bgm  = like.obs().npred()
			
			# Assess quality based on a comparison between Npred and Nevents
			quality_bgm = npred_bgm-nevents

			# Write background fit results
			if self.logExplicit():
				self.log.header3("Background model fit")
				self.log.parformat("log likelihood")
				self.log(LogL_bgm)
				self.log("\n")
				self.log.parformat("Number of predicted events")
				self.log(npred_bgm)
				self.log("\n")
				self.log.parformat("Fit quality")
				self.log(quality_bgm)
				self.log("\n")

			# Start over if the fit quality was bad
			if abs(quality_bgm) > 3.0:
				if self.logExplicit():
					self.log("Fit quality outside required range. Start over.\n")
				continue
			
			# Write model fit results
			if self.logExplicit():
				for model in result_bgm:
					self.log.parformat("Model")
					self.log(model.name())
					self.log("\n")
					for par in model:
						self.log(str(par)+"\n")
			
			# Fit background and test source
			sim.models(full_model)
			like       = obsutils.fit(sim, log=self.m_log, debug=self.m_debug)
			result_all = like.obs().models().copy()
			LogL_all   = like.opt().value()
			npred_all  = like.obs().npred()
			ts         = 2.0*(LogL_bgm-LogL_all)

			# Assess quality based on a comparison between Npred and Nevents
			quality_all = npred_all-nevents

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
				self.log.parformat("Fit quality")
				self.log(quality_all)
				self.log("\n")
				#
				for model in result_all:
					self.log.parformat("Model")
					self.log(model.name())
					self.log("\n")
					for par in model:
						self.log(str(par)+"\n")

			# Start over if the fit quality was bad
			if abs(quality_all) > 3.0:
				if self.logExplicit():
					self.log("Fit quality outside required range. Start over.\n")
				continue

			# Start over if TS was non-positive
			if ts <= 0.0:
				if self.logExplicit():
					self.log("Non positive TS. Start over.\n")
				continue

			# Get fitted Crab, photon and energy fluxes
			cflux     = result_all["Test"]['Prefactor'].value()
			cflux_err = result_all["Test"]['Prefactor'].error()
			pflux     = cast_GModelSky(result_all["Test"]).spectral().flux(emin, emax)
			eflux     = cast_GModelSky(result_all["Test"]).spectral().eflux(emin, emax)

			# Compute differential sensitivity
			diffSens = pflux / erg_width * erg_mean*erg_mean

			# Compute flux correction factor based on average TS
			correct = 1.0
			if self.m_use_ts:
				if ts > 0:
					correct = math.sqrt(self.m_ts_thres/ts)
			else:
				if cflux > 0 and cflux_err > 0:
					correct = math.sqrt(self.m_ts_thres)/(cflux/cflux_err)
			
			# Compute extrapolated fluxes
			flux     = correct * cflux
			pflux    = correct * pflux
			eflux    = correct * eflux
			diffSens = correct * diffSens
			flux_value.append(flux)
			pflux_value.append(pflux)
			eflux_value.append(eflux)
			diffSens_value.append(diffSens)
			
			# Write background and test source fit results
			if self.logExplicit():
				self.log.parformat("Photon flux")
				self.log(pflux)
				self.log(" ph/cm2/s\n")
				self.log.parformat("Energy flux")
				self.log(eflux)
				self.log(" erg/cm2/s\n")
				self.log.parformat("Crab flux")
				self.log(cflux*1000.0)
				self.log(" mCrab\n")
				self.log.parformat("Differential sensitivity")
				self.log(diffSens)
				self.log(" erg/cm2/s\n")
				for model in result_all:
					self.log.parformat("Model")
					self.log(model.name())
					self.log("\n")
					for par in model:
						self.log(str(par)+"\n")
			elif self.logTerse():
				self.log.parformat("Iteration "+str(iter))
				self.log("TS=")
				self.log(ts)
				self.log(" ")
				self.log("corr=")
				self.log(correct)
				self.log("  ")
				self.log(pflux)
				self.log(" ph/cm2/s = ")
				self.log(eflux)
				self.log(" erg/cm2/s = ")
				self.log(cflux*1000.0)
				self.log(" mCrab = ")
				self.log(diffSens)
				self.log(" erg/cm2/s\n")
			
			# Compute sliding average of extrapolated fitted prefactor,
			# photon and energy flux. This damps out fluctuations and
			# improves convergence
			flux     = 0.0
			pflux    = 0.0
			eflux    = 0.0
			diffSens = 0.0
			num      = 0.0
			for k in range(self.m_num_avg):
				inx = len(flux_value) - k - 1
				if inx >= 0:
					flux     += flux_value[inx]
					pflux    += pflux_value[inx]
					eflux    += eflux_value[inx]
					diffSens += diffSens_value[inx]
					num      += 1.0
			flux     /= num
			pflux    /= num
			eflux    /= num
			diffSens /= num
			
			# Compare average flux to last average
			if iter > self.m_num_avg:
				if test_flux > 0:
					ratio = flux/test_flux
					
					# We have 2 convergence criteria:
					# 1. The average flux does not change
					# 2. The flux correction factor is small
					if ratio   >= 0.99 and ratio   <= 1.01 and \
					   correct >= 0.9  and correct <= 1.1:
						if self.logTerse():
							self.log(" Converged ("+str(ratio)+")\n")
						break
				else:
					if self.logTerse():
						self.log(" Flux is zero.\n")
					break
			
			# Use average for next iteration
			test_flux = flux
			
			# Exit loop if number of trials exhausted
			if (iter >= self.m_max_iter):
				if self.logTerse():
					self.log(" Test ended after "+str(self.m_max_iter)+" iterations.\n")
				break

		# Write fit results
		if self.logTerse():
			self.log.header3("Fit results")
			self.log.parformat("Test statistics")
			self.log(ts)
			self.log("\n")
			self.log.parformat("Photon flux")
			self.log(pflux)
			self.log(" ph/cm2/s\n")
			self.log.parformat("Energy flux")
			self.log(eflux)
			self.log(" erg/cm2/s\n")
			self.log.parformat("Crab flux")
			self.log(cflux*1000.0)
			self.log(" mCrab\n")
			self.log.parformat("Differential sensitivity")
			self.log(diffSens)
			self.log(" erg/cm2/s\n")
			self.log.parformat("Number of simulated events")
			self.log(nevents)
			self.log("\n")
			self.log.header3("Background and test source model fitting")
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
			
		# Store result
		result = {'loge': loge, 'emin': emin.TeV(), 'emax': emax.TeV(), \
		          'crab': flux, 'pflux': pflux, 'eflux': eflux, \
				  'diffSens': diffSens}
		
		# Return result
		return result


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
	"""
	CTA ctlike-based sensitivity estimator.
	"""
	# Create instance of application
	app = cssens(sys.argv)
	
	# Run application
	app.execute()
	