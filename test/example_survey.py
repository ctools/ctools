#! /usr/bin/env python
# ===========================================================================================#
# This script simulates CTA surveys.
#
# If matplotlib is installed, the spectrum will be displayed on the screen. 
# ===========================================================================================#
from ctatools import *
from gammalib import *
from math import *
import os
import glob
import sys


# =========== #
# Plot counts #
# =========== #
def plot_counts(obs):
	"""
	Plot counts.
	"""
	# Only proceed if matplotlib is available
	if 1 == 1:
	#try:
		# Import matplotlib
		import matplotlib.pyplot as plt
		
		# Set legend fontsize
		params = {'legend.fontsize': 10}
		plt.rcParams.update(params)
		
		# Set plot styles
		styles = ['b-', 'g-', 'y-', 'n-']
		
		# Dump header
		print ""
		print "Make plots (using matplotlib):"
		print "=============================="
		
		# Create figure 1
		plt.figure(1,figsize=(12,6))
		plt.subplots_adjust(hspace=.7)
		
		# Create subplot 1
		plt.subplot(121)
		plt.title("Spectrum (summed over all pixels)")
		
		# Loop over observations
		for run in obs:
		
			# Get event list
			list = cast_GCTAEventList(run.events())
			
			# Create energy axis
			ebounds = GEbounds()
			emin    = GEnergy()
			emax    = GEnergy()
			emin.TeV(0.1)
			emax.TeV(100.0)
			ebounds.setlog(emin, emax, 10)
			energy = [ebounds.elogmean(i).TeV() for i in range(ebounds.size())]
			
			# Create spectrum
			print "Extract data:"
			counts = [0.0 for i in range(ebounds.size())]
			print list.size()
			for atom in list:
				index         = ebounds.index(atom.energy())
				counts[index] = counts[index] + 1.0
			
			# Create error bars
			error = [sqrt(c) for c in counts]
			
			# Plot spectrum
			plt.loglog(energy, counts, 'ro', label='data')
			#plt.errorbar(energy, counts, error, fmt=None, ecolor='r')
			
			# Extract models
			print "Extract models:"
			sum_model = [0.0 for i in range(ebounds.size())]
			for k, m in enumerate(obs.models()):
				print "- "+m.name()
				model = [0.0 for i in range(ebounds.size())]
				for atom in list:
					index        = ebounds.index(atom.energy())
					prob         = m.eval(atom, run)
					model[index] = model[index] + prob * atom.size()
				for i in range(ebounds.size()):
					sum_model[i] = sum_model[i] + model[i]
				#plt.loglog(energy, model, styles[k], label=m.name())
			#plt.loglog(energy, sum_model, 'r-', label='total')
		
		# Put labels
		plt.xlabel("Energy (TeV)")
		plt.ylabel("Counts")
		plt.legend(loc="lower left")
		
		# Create subplot 2
		plt.subplot(122)
		plt.title("Offset (summed over all energies)")
		
		# Set Crab direction
		crab = GSkyDir()
		crab.radec_deg(83.63, 22.01)
		
		# Loop over observations
		for run in obs:
		
			# Get event cube
			list = cast_GCTAEventList(run.events())
			
			# Create offset histogram
			print "Extract data:"
			nx       = 30
			doffset2 = 0.01
			offset2  = [(i+0.5)*doffset2 for i in range(nx)]
			counts   = [0.0  for i in range(nx)]
			for atom in list:
				off   = atom.dir().dist_deg(crab)
				off2  = off*off
				index = int(off2/doffset2)
				if index < nx:
					counts[index] = counts[index] + 1.0
			
			# Create error bars
			error = [sqrt(c) for c in counts]
			
			# Plot distribution
			#plt.semilogy(offset2, counts, 'ro', label='data')
			
			# Extract models
			print "Extract models:"
			sum_model = [0.0 for i in range(nx)]
			for k, m in enumerate(obs.models()):
				print "- "+m.name()
				model = [0.0 for i in range(nx)]
				for atom in list:
					off   = atom.dir().dist_deg(crab)
					off2  = off*off
					index = int(off2/doffset2)
					if index < nx:
						prob         = m.eval(atom, run)
						model[index] = model[index] + prob * atom.size()
				for i in range(nx):
					sum_model[i] = sum_model[i] + model[i]
				#plt.plot(offset2, model, styles[k], label=m.name())
			#plt.plot(offset2, sum_model, 'r-', label='total')
			#plt.ylim(ymin=0.1)
		
		# Put labels
		plt.xlabel("Offset (deg^2)")
		plt.ylabel("Counts")
		#plt.legend(loc="upper right")
		
		# Show counts spectra
		plt.show()
		
	#except:
	#	print "Matplotlib is not (correctly) installed on your system. No counts spectra are shown."

	# Return
	return


# ===================== #
# Simulate observations #
# ===================== #
def sim_obs(obs):
	"""
	Simulate events for all observations in the container.
	"""
	# Simulate events
	sim = ctobssim(obs)
	sim.logFileOpen()  # We need this to explicitely open the log file in Python mode
	sim.run()
	
	# Extract copy of observations (we need a copy here as sim goes out of
	# scope once we leave the function)
	obs = sim.obs().copy()
	
	# Return observations
	return obs


# ================ #
# Fit observations #
# ================ #
def fit_obs(obs):
	"""
	Perform maximum likelihood fitting of observations.
	"""
	# Perform maximum likelihood fitting
	like = ctlike(obs)
	like.logFileOpen()  # We need this to explicitely open the log file in Python mode
	like.run()
	
	# Extract copy of observations (we need a copy here as like goes out of
	# scope once we leave the function)
	#obs = like.obs().copy()
	
	#print like.opt()
	
	# Return observations
	return like


# ================= #
# Create counts map #
# ================= #
def make_counts_map(obs):
	"""
	Make counts map by combining the events of all observations.
	"""
	# Allocate counts map
	nxpix  = 600
	nypix  = 200
	cntmap = GSkymap("CAR", "GAL", 0.0, 0.0, -0.05, 0.05, nxpix, nypix, 1)
	
	# Fill all observations
	for run in obs:
		
		# Loop over all events
		for event in run.events():
			
			# Cast to CTA event
			event = cast_GCTAEventAtom(event)
			
			# Determine sky pixel
			skydir = event.dir().skydir()
			pixel  = cntmap.dir2pix(skydir)
			
			# Set pixel
			cntmap[pixel] += 1.0
	
	# Save sky map
	cntmap.save("survey_cntmap.fits", True)
	
	# Return counts map
	return cntmap


# ================ #
# Create model map #
# ================ #
def make_model_map(obs):
	"""
	Make model map by combining all observations.
	"""
	# Allocate model map
	nxpix  = 600
	nypix  = 200
	modmap = GSkymap("CAR", "GAL", 0.0, 0.0, -0.05, 0.05, nxpix, nypix, 1)
	
	# Set energy and time
	energy  = GEnergy()
	time    = GTime()
	instdir = GCTAInstDir()
	energy.TeV(0.1)
	
	# Loop over all map pixels
	for pixel in range(modmap.npix()):
		
		# Get sky direction
		skydir = modmap.pix2dir(pixel)
		instdir.skydir(skydir)
		
		# Create event atom for map pixel
		atom = GCTAEventAtom()
		atom.dir(instdir)
		atom.energy(energy)
		atom.time(time)
		
		# Initialise model value
		value = 0.0
		
		# Loop over all observations
		for run in obs:
			value += obs.models().eval(atom, run)
		
		# Set map value
		modmap[pixel] = value
	
	# Save sky map
	modmap.save("survey_modmap.fits", True)
	
	# Return model map
	return modmap


# ======================= #
# Set one CTA observation #
# ======================= #
def set_one_obs(pntdir, tstart=0.0, duration=1800.0, emin=0.1, emax=100.0, rad=3.0, \
                irf="kb_E_50h_v3", caldb="irf"):
	"""
	Parameters:
	 pntdir   - Pointing direction
	Keywords:
	 tstart   - Start time (seconds)
	 duration - Duration of observation (seconds)
	 emin     - Minimum event energy (TeV)
	 emax     - Maximum event energy (TeV)
	 rad      - ROI radius (deg)
	 irf      - Instrument response function
	 caldb    - Calibration database path
	"""
	# Allocate CTA observation
	obs = GCTAObservation()
	
	# Set pointing direction
	pnt    = GCTAPointing()
	pnt.dir(pntdir)
	obs.pointing(pnt)
	
	# Set ROI
	roi     = GCTARoi()
	instdir = GCTAInstDir()
	instdir.skydir(pntdir)
	roi.centre(instdir)
	roi.radius(rad)
	
	# Set GTI
	gti   = GGti()
	start = GTime()
	stop  = GTime()
	start.met(tstart)
	stop.met(tstart+duration)
	gti.append(start, stop)
	
	# Set energy boundaries
	ebounds = GEbounds()
	e_min   = GEnergy()
	e_max   = GEnergy()
	e_min.TeV(emin)
	e_max.TeV(emax)
	ebounds.append(e_min, e_max)

	# Allocate event list
	events = GCTAEventList()
	events.roi(roi)
	events.gti(gti)
	events.ebounds(ebounds)
	obs.events(events)
	
	# Set instrument response
	obs.response(irf, caldb)
	
	# Return observation
	return obs


# ======================================== #
# Add CTA background model to observations #
# ======================================== #
def add_background_model(obs):
	"""
	Add standard CTA background model to observations container.
	
	We use a simple power law here, scaled to Konrad's E configuration
	performance table. The model needs still to be validated.
	"""
	# Recover models from observation
	models = obs.models()
	
	# Define background model
	bgd_radial   = GCTAModelRadialGauss(3.0)
	bgd_spectrum = GModelSpectralPlaw(61.8, -1.85)
	bgd_spectrum["Prefactor"].scale(1.0e-6)
	bgd_spectrum["PivotEnergy"].value(1.0)
	bgd_spectrum["PivotEnergy"].scale(1.0e6)
	bgd_model = GCTAModelRadialAcceptance(bgd_radial, bgd_spectrum)
	bgd_model.name("Background")
	bgd_model.instruments("CTA")
	
	# Add background model to container
	models.append(bgd_model)
	
	# Put container back in observation container
	obs.models(models)
	
	# Return observation container
	return obs


# ================= #
# Set Crab spectrum #
# ================= #
def crab_spec():
	"""
	Set Crab spectrum based on MAGIC observations
	(Albert et al. 2008, ApJ, 674, 1037)
	"""
	# Set parameters
	spectrum = GModelSpectralPlaw(5.7, -2.48)
	spectrum["Prefactor"].scale(1.0e-16)
	spectrum["PivotEnergy"].value(0.3)
	spectrum["PivotEnergy"].scale(1.0e6)
	
	# Return spectrum
	return spectrum


# =============================== #
# Setup single observation survey #
# =============================== #
def survey_single():
	"""
	Creates a single observation survey for test purposes.
	"""
	# Allocate observation container
	obs = GObservations()
	
	# Set single pointing at galactic centre
	pntdir = GSkyDir()
	pntdir.lb_deg(0.0, 0.0)
	run = set_one_obs(pntdir)
	obs.append(run)
	
	# Define single point source with Crab flux at galactic centre
	center = GSkyDir()
	center.lb_deg(0.0, 0.0)
	point_spatial  = GModelSpatialPtsrc(center)
	point_spectrum = crab_spec()
	point          = GModelPointSource(point_spatial, point_spectrum)
	point.name('GC source')
	
	# Create model container
	models = GModels()
	models.append(point)
	obs.models(models)
	
	# Return observation container
	return obs


# =========================== #
# Setup Galactic plane survey #
# =========================== #
def survey_gplane(lrange=10, lstep=2):
	"""
	Creates a single observation survey for test purposes.
	
	Keywords:
	 lrange - Longitude range (integer deg)
	 lstep  - Longitude step size (integer deg)
	"""
	# Allocate observation container
	obs = GObservations()
	
	# Loop over longitudes
	for l in range(-lrange,lrange+lstep,lstep):
		
		# Set pointing
		pntdir = GSkyDir()
		pntdir.lb_deg(l, 0.0)
		run = set_one_obs(pntdir)
		obs.append(run)
	
	# Define single point source with Crab flux at galactic centre
	center = GSkyDir()
	center.lb_deg(0.0, 0.0)
	point_spatial  = GModelSpatialPtsrc(center)
	point_spectrum = crab_spec()
	point          = GModelPointSource(point_spatial, point_spectrum)
	point.name('GC source')
	
	# Create model container
	models = GModels()
	models.append(point)
	obs.models(models)
	
	# Return observation container
	return obs


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	CTA survey simulation and analysis script.
	"""
	# Initialise flags
	need_help = False
	
	# Test for command line arguments
	print sys.argv[0]
	if (len(sys.argv) > 1):
		if sys.argv[1] == "-h":
			need_help = True
		else:
			need_help = True

	# Print help if needed and exit
	if need_help:
		print "Usage: example_survey.py [OPTIONS]"
		print "     -h       Display this usage message"
		sys.exit()
	
	# Dump header
	print "**********************"
	print "* Simulate CTA survey *"
	print "**********************"
	
	# Remove any existing result files
	list = [glob.glob("*.fits"), glob.glob("*.log"), glob.glob("*.xml")]
	for files in list:
		for file in files:
			os.remove(file)
	
	# Setup single observation survey
	obs = survey_single()
	obs = survey_gplane()
	
	# Add background model
	obs = add_background_model(obs)
	
	# Simulate events
	print "Simulate events"
	obs = sim_obs(obs)
	
	# Make counts map
	print "Make counts map"
	cntmap = make_counts_map(obs)
	
	# Fit observations
	print "Fit observatins"
	like = fit_obs(obs)
	#print like.opt()
	#print like.obs().models()
	
	# Make model map
	print "Make model map (this step will take some time)"
	modmap = make_model_map(obs)
	
	# Show fit results
	#plot_counts(like.obs())
