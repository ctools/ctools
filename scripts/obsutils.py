# ==========================================================================
# This script provides a number of functions that are useful for handling
# CTA data.
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
from math import *
#import os
#import glob
#import sys
#import csv
#import math
#import random


# ===================== #
# Simulate observations #
# ===================== #
def sim_obs(obs, log=False, seed=0, nbins=0, binsz=0.05, npix=200):
	"""
	Simulate events for all observations in the container.
	
	Parameters:
	 obs   - Observation container
	Keywords:
	 log   - Create log file(s)
	 seed  - Seed value for simulations (default: 0)
	 nbins - Number of energy bins (default: 0=unbinned)
	 binsz - Pixel size for binned simulation (deg/pixel)
	 npix  - Number of pixels in X and Y for binned simulation
	"""
	# Allocate ctobssim application and set parameters
	sim = ctobssim(obs)
	sim['seed'].integer(seed)
		
	# Optionally open the log file
	if log:
		sim.logFileOpen()
	
	# Run ctobssim application. This will loop over all observations in the
	# container and simulation the events for each observation. Note that
	# events are not added together, they still apply to each observation
	# separately.
	sim.run()
	
	# Binned option?
	if nbins > 0:
	
		# Determine common energy boundaries for observations
		emin = None
		emax = None
		for run in sim.obs():
			run_emin = run.events().ebounds().emin().TeV()
			run_emax = run.events().ebounds().emax().TeV()
			if emin == None:
				emin = run_emin
			elif run_emin > emin:
				emin = run_emin
			if emax == None:
				emax = run_emax
			elif run_emax > emax:
				emax = run_emax
	
		# Allocate ctbin application and set parameters
		bin = ctbin(sim.obs())
		bin["emin"].real(emin)
		bin["emax"].real(emax)
		bin["enumbins"].integer(nbins)
		bin["nxpix"].integer(npix)
		bin["nypix"].integer(npix)
		bin["binsz"].real(binsz)
		bin["coordsys"].string("GAL")
		bin["proj"].string("TAN")
		
		# Optionally open the log file
		if log:
			bin.logFileOpen()

		# Run ctbin application. This will loop over all observations in
		# the container and bin the events in counts maps
		bin.run()
		
		# Make a deep copy of the observation that will be returned
		# (the ctbin object will go out of scope one the function is
		# left)
		obs = bin.obs().copy()
		
	else:

		# Make a deep copy of the observation that will be returned
		# (the ctobssim object will go out of scope one the function is
		# left)
		obs = sim.obs().copy()
	
	# Return observation container
	return obs


# ================ #
# Fit observations #
# ================ #
def fit_obs(obs, log=False):
	"""
	Perform maximum likelihood fitting of observations in the container.
	
	Parameters:
	 obs - Observation container
	Keywords:
	 log - Create log file(s)
	"""
	# Allocate ctlike application
	like = ctlike(obs)
	
	# Optionally open the log file
	if log:
		like.logFileOpen()
		like["debug"].boolean(True)

	# Run ctlike application.
	like.run()
	
	# Return observations
	return like


# ================= #
# Create counts map #
# ================= #
def make_counts_map(obs, proj="TAN", coord="GAL", xval=0.0, yval=0.0, \
                         nxpix=600, nypix=200, pixsize=0.05, \
						 outname="survey_cntmap.fits"):
	"""
	Creates a counts map by combining the events of all observations.
	The counts map will be a summed map over all energies.
	
	Parameters:
	 obs - Observation container
	Keywords:
	 xval    - Reference longitude value [deg]
	 yval    - Reference latitude value [deg]
	 nxpix   - Number of pixels in longitude axis
	 nypix   - Number of pixels in latitude axis
	 pixsize - Pixel size [deg]
	 outname - Counts map FITS filename
	"""
	# Allocate counts map
	cntmap = GSkymap(proj, coord, xval, yval, -pixsize, pixsize, nxpix, nypix, 1)
	
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
	
	# Save sky map. The clobber flag is set to True, so any existing FITS
	# file will be overwritten.
	cntmap.save(outname, True)
	
	# Return counts map
	return cntmap


# ================ #
# Create model map #
# ================ #
def make_model_map(obs, eref=0.1, proj="TAN", coord="GAL", xval=0.0, yval=0.0, \
						nxpix=600, nypix=200, pixsize=0.05, \
						outname="survey_modmap.fits"):
	"""
	Make model map for a given reference energy by combining all observations.
	The model map will be evaluated for a given reference energy 'eref' and will
	be given in units of [counts/(sr MeV s)].
	
	Parameters:
	 obs - Observation container
	Keywords:
	 eref    - Reference energy for which model is created [TeV]
	 xval    - Reference longitude value [deg]
	 yval    - Reference latitude value [deg]
	 nxpix   - Number of pixels in longitude axis
	 nypix   - Number of pixels in latitude axis
	 pixsize - Pixel size [deg]
	 outname - Counts map FITS filename
	"""
	# Allocate model map
	modmap = GSkymap(proj, coord, xval, yval, -pixsize, pixsize, nxpix, nypix, 1)
	
	# Set reference energy, time and direction. The time is not initialised and is
	# in fact not used (as the IRF is assumed to be time independent for now).
	# The sky direction is set later using the pixel values.
	energy  = GEnergy()
	time    = GTime()
	instdir = GCTAInstDir()
	energy.TeV(eref)
	
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
	modmap.save(outname, True)
	
	# Return model map
	return modmap


# ======================= #
# Set one CTA observation #
# ======================= #
def set_one_obs(pntdir, tstart=0.0, duration=1800.0, \
                emin=0.1, emax=100.0, rad=5.0, \
                irf="kb_E_50h_v3", caldb="/usr/local/gamma/share/caldb/cta"):
	"""
	Returns a single CTA observation. By looping over this function we can
	add CTA observations to the observation container.
	
	Parameters:
	 pntdir   - Pointing direction
	Keywords:
	 tstart   - Start time [seconds]
	 duration - Duration of observation [seconds]
	 emin     - Minimum event energy [TeV]
	 emax     - Maximum event energy [TeV]
	 rad      - ROI radius used for analysis [deg]
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
	#print obs
	
	# Return observation
	return obs


# ================== #
# Setup source model #
# ================== #
def add_source_model(models, source):
	"""
	This function adds one source model to the observations container.

	Parameters:
	 obs    - Observations container
	 source - Source to be added to container
	"""
	# Append source to model container
	models.append(source)
	
	# Return model container
	return models


# ======================================== #
# Add CTA background model to observations #
# ======================================== #
def add_background_model(models):
	"""
	Add standard CTA background model to observations container.
	We use a simple power law here, scaled to Konrad's E configuration
	performance table.

	Parameters:
	 obs - Observation container
	"""
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
	
	# Return model container
	return models


# ==================== #
# Set Crab-like source #
# ==================== #
def crab_like(name, l, b, flux=1.0, type="point", sigma=1.0, radius=1.0, width=0.1, \
              fitpos=True):
	"""
	Returns a single source with Crab-like spectrum. The source flux
	can be scaled in Crab units.

	Parameters:
	 name   - Unique source name
	 l      - Galactic longitude of source location [deg]
	 b      - Galactic latitude of source location [deg]
	Keywords:
	 flux   - Source flux in Crab units
	 type   - Source type ("point", "gauss", "disk", "shell")
	 sigma  - Gaussian sigma (for type="gauss")
	 radius - Disk or shell inner radius [deg] (for type="disk" and type="shell")
	 width  - Shell width [deg] (for type="shell")
	 fitpos - Specifies if position should be fitted or not
	"""
	# Set source location
	location = GSkyDir()
	location.lb_deg(l, b)
	
	# Set source spectrum
	spectrum = crab_spec(flux=flux)

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
		print "ERROR: Unknown source type '"+type+"'."
		return None
	
	# Set source name
	source.name(name)
	
	# Return source
	return source


# ================= #
# Set Crab spectrum #
# ================= #
def crab_spec(flux=1.0):
	"""
	Set Crab spectrum based on MAGIC observations
	(Albert et al. 2008, ApJ, 674, 1037)

	Keywords:
	 flux - Source flux in Crab units
	"""
	# Set parameters
	spectrum = GModelSpectralPlaw(flux, -2.48)
	spectrum["Prefactor"].scale(5.7e-16)
	spectrum["PivotEnergy"].value(0.3)
	spectrum["PivotEnergy"].scale(1.0e6)
	
	# Return spectrum
	return spectrum


# ========================= #
# Set any point-like source #
# ========================= #
def any_pointlike(name, l, b, flux=1.0, index=-2.48, fitpos=True, fitidx=True):
	"""
	Returns a single point-source with power-law spectrum.
	The source flux can be scaled in Crab units.

	Parameters:
	 name   - Unique source name
	 l      - Galactic longitude of source location [deg]
	 b      - Galactic latitude of source location [deg]
	Keywords:
	 flux   - Source flux in Crab units
	 fitpos - Specifies if position should be fitted or not
	 fitidx - Specifies if index should be fitted or not
	"""
	# Set source location
	location = GSkyDir()
	location.lb_deg(l, b)
	
	
	# Set source spatial model
	spatial = GModelSpatialPtsrc(location)
	
	# Set source spectral model
	spectrum = GModelSpectralPlaw(flux, index)
	spectrum["Prefactor"].scale(5.7e-16)
	spectrum["PivotEnergy"].value(0.3)
	spectrum["PivotEnergy"].scale(1.0e6)
	
	# Set fit properties
	if fitpos:
		spatial[0].free()
		spatial[1].free()
	else:
		spatial[0].fix()
		spatial[1].fix()
	
	if fitidx:
		spectrum["Index"].free()
	else:
		spectrum["Index"].fix()	

	# Assign spatial and spectral models to source object
	source  = GModelPointSource(spatial, spectrum)
	
	# Set source name
	source.name(name)
	
	# Return source
	return source


# =============================== #
# Setup single observation survey #
# =============================== #
def survey_single(caldb, irf, lpnt=0.0, bpnt=0.0):
	"""
	This is an example of how to create a pointing strategy. Here, we
	only have a single CTA pointing, so this is a very basic pointing
	strategy.

	Parameters:
	 caldb - Calibration database
	 irf   - Response function
	 lpnt  - Galactic longitude of pointing [deg]
	 bpnt  - Galactic latitude of pointing [deg]
	"""
	# Allocate observation container
	obs = GObservations()
	
	# Set single pointing at galactic centre
	pntdir = GSkyDir()
	pntdir.lb_deg(lpnt, bpnt)
	run = set_one_obs(pntdir, caldb=caldb, irf=irf)
	obs.append(run)
	
	# Return observation container
	return obs


# =========================== #
# Setup Galactic plane survey #
# =========================== #
def survey_gplane(caldb, irf, bpnt=0.0, lrange=10, lstep=2):
	"""
	Creates a single observation survey for test purposes.
	
	Keywords:
	 bpnt   - Latitude value for scan [deg]
	 lrange - Longitude range (integer deg)
	 lstep  - Longitude step size (integer deg)
	"""
	# Allocate observation container
	obs = GObservations()
	
	# Loop over longitudes
	for l in range(-lrange,lrange+lstep,lstep):
		
		# Set pointing
		pntdir = GSkyDir()
		pntdir.lb_deg(l, bpnt)
		run = set_one_obs(pntdir, caldb=caldb, irf=irf)
		obs.append(run)
	
	# Return observation container
	return obs


# =========================== #
# Setup Galactic plane survey #
# =========================== #
def survey_square(caldb, irf, lrange=2.0, lstep=2.0, duration=1800.0, \
                  emin=0.1, emax=100.0):
	"""
	Creates a single observation survey over a square grid.
	For the moment, the grid is centred on (l,b)=(0,0)
	
	Arguments:
	 lrange - Longitude/latitude range (float deg)
	 lstep  - Longitude/latitude step size (float deg)
	 duration/emin/emax - See set_one_obs function
	"""
	# Allocate observation container
	obs = GObservations()
	
	# Compute number of points along each axis
	num_pnt= int(round(2.0*lrange/lstep)+1.0)
	
	# Loop over longitudes and latitudes
	for i in range(0,num_pnt):
		for j in range(0,num_pnt):
		
			# Set pointing
			l= -lrange+i*lstep
			b= -lrange+j*lstep
			pntdir = GSkyDir()
			pntdir.lb_deg(l, b)
			run = set_one_obs(pntdir, caldb=caldb, irf=irf, duration=duration, \
			                  emin=emin, emax=emax)
			obs.append(run)
	
	# Return observation container
	return obs


# =========================== #
# Setup Galactic plane survey #
# =========================== #
def survey_galplane(caldb, irf, lmin=-2.0, lstep=2.0, nstep=3, duration=1800.0, \
                    emin=0.1, emax=100.0):
	"""
	Creates a single observation survey along the galactic plane
	
	Arguments:
	 lmin   - Starting longitude(float deg)
	 lstep  - Longitude step size (float deg)
	 nstep  - Number of steps/pointings
	 duration/emin/emax - See set_one_obs function
	"""
	# Allocate observation container
	obs = GObservations()
		
	# Loop over longitudes
	for i in range(0,nstep):
		
		# Set pointing
		l= lmin+i*lstep
		b= 0.0
		pntdir = GSkyDir()
		pntdir.lb_deg(l, b)
		run = set_one_obs(pntdir, caldb=caldb, irf=irf, duration=duration, \
		                  emin=emin, emax=emax)
		obs.append(run)
	
	# Return observation container
	return obs

# ===================================== #
# Setup set of point source Crab models #
# ===================================== #
def add_point_sources(obs):
	"""
	Add point-source model made of Crab-like sources.
	"""
	# Add point-source Crabs
	obs = add_source_model(obs, crab_like("Crab", 0.0, 0.0))
	obs = add_source_model(obs, crab_like("Crab 1deg off", 1.0, 0.0))
	obs = add_source_model(obs, crab_like("Crab 2deg off", 2.0, 0.0))
	obs = add_source_model(obs, crab_like("Crab 3deg off", 3.0, 0.0))
	obs = add_source_model(obs, crab_like("Crab 4deg off", 4.0, 0.0))
	obs = add_source_model(obs, crab_like("Crab 5deg off", 5.0, 0.0))
	
	# Return
	return obs


# ================= #
# Setup disk models #
# ================= #
def add_disk_sources(obs):
	"""
	Add disk-source model made of Crab-like sources.
	"""
	# Add disk Crabs
	obs = add_source_model(obs, crab_like("Disk 0deg off",  0.0, 0.0, type="disk", radius=0.3))
	obs = add_source_model(obs, crab_like("Disk 1deg off", -1.0, 0.0, type="disk", radius=0.3))
	obs = add_source_model(obs, crab_like("Disk 2deg off", -2.0, 0.0, type="disk", radius=0.3))
	obs = add_source_model(obs, crab_like("Disk 3deg off", -3.0, 0.0, type="disk", radius=0.3))
	#obs = add_source_model(obs, crab_like("Disk 4deg off", -4.0, 0.0, type="disk", radius=0.3))
	#obs = add_source_model(obs, crab_like("Disk 5deg off", -5.0, 0.0, type="disk", radius=0.3))
	
	# Return
	return obs


# ================== #
# Setup Gauss models #
# ================== #
def add_gauss_sources(obs):
	"""
	Add Gauss-source model made of Crab-like sources.
	"""
	# Add Gauss Crabs
	obs = add_source_model(obs, crab_like("Gauss 0deg off",  0.0, 0.0, type="gauss", sigma=0.3))
	obs = add_source_model(obs, crab_like("Gauss 1deg off", -1.0, 0.0, type="gauss", sigma=0.3))
	obs = add_source_model(obs, crab_like("Gauss 2deg off", -2.0, 0.0, type="gauss", sigma=0.3))
	obs = add_source_model(obs, crab_like("Gauss 3deg off", -3.0, 0.0, type="gauss", sigma=0.3))
	#obs = add_source_model(obs, crab_like("Gauss 4deg off", -4.0, 0.0, type="gauss", sigma=0.3))
	#obs = add_source_model(obs, crab_like("Gauss 5deg off", -5.0, 0.0, type="gauss", sigma=0.3))

	# Return
	return obs


# ================== #
# Setup Shell models #
# ================== #
def add_shell_sources(obs):
	"""
	Add shell-source model made of Crab-like sources.
	"""
	# Add Shell Crabs
	obs = add_source_model(obs, crab_like("Shell 0deg off",  0.0, 0.0, type="shell", radius=0.3, width=0.1))
	obs = add_source_model(obs, crab_like("Shell 1deg off", -1.0, 0.0, type="shell", radius=0.3, width=0.1))
	obs = add_source_model(obs, crab_like("Shell 2deg off", -2.0, 0.0, type="shell", radius=0.3, width=0.1))
	#obs = add_source_model(obs, crab_like("Shell 3deg off", -3.0, 0.0, type="shell", radius=0.3, width=0.1))
	#obs = add_source_model(obs, crab_like("Shell 4deg off", -4.0, 0.0, type="shell", radius=0.3, width=0.1))
	#obs = add_source_model(obs, crab_like("Shell 5deg off", -5.0, 0.0, type="shell", radius=0.3, width=0.1))
	
	# Return
	return obs


# ==================== #
# Setup big disk model #
# ==================== #
def add_big_disk_source(obs):
	"""
	Add big disk-source model.
	"""
	# Add one big offset disk Crab
	obs = add_source_model(obs, crab_like("Big Disk", 1.0, 1.0, type="disk", radius=1.0))
	
	# Return
	return obs


# ===================== #
# Setup big Gauss model #
# ===================== #
def add_big_gauss_source(obs):
	"""
	Add big gauss-source model.
	"""
	# Add one big offset Gauss Crab
	obs = add_source_model(obs, crab_like("Big Gauss", -1.0, 1.0, type="gauss", sigma=1.0))
	
	# Return
	return obs


# ===================== #
# Setup big shell model #
# ===================== #
def add_big_shell_source(obs):
	"""
	Add big shell-source model.
	"""
	# Add one big offset Shell Crab
	obs = add_source_model(obs, crab_like("Big Shell", -1.0, -1.0, type="shell", radius=1.0, width=0.1))
	
	# Return
	return obs


# ================================================= #
# Setup set of point source for the survey analysis #
# ================================================= #
def add_survey_sources(obs):
	"""
	Add point-source model made of Crab-like sources.
	"""
	# Add point-source Crabs
	obs = add_source_model(obs, crab_like("Source 1", 0.0, 0.0, flux=0.1, type="point"))
	obs = add_source_model(obs, crab_like("Source 2", 1.0, 0.0, flux=0.2, type="point"))
	obs = add_source_model(obs, crab_like("Source 3", -1.0, 0.0, flux=0.5, type="point"))
	obs = add_source_model(obs, crab_like("Source 4", 0.0, 1.0, flux=1.0, type="point"))
	obs = add_source_model(obs, crab_like("Source 5", 0.0, -1.0, flux=2.0, type="point"))
	obs = add_source_model(obs, crab_like("Source 6", 0.0, 4.0, flux=5.0, type="point"))
	
	# Return
	return obs


# ================================================= #
# Setup set of point source for the survey analysis #
# ================================================= #
def add_random_sources(model):
	"""
	Add point-sources uniformly distributed over the sky, with a parameterised distribution in flux.
	
	Arguments
	- model: a source model container to hold the source model over the loop
	
	Notes
	- Contrary to other add_source functions, this one works on a model container and not on an observation container
	- This is because we need to reuse the same random source model in a loop over pointing strategies 
	"""
	# Define source population characteristics
	flux_min=0.001
	flux_max=10.0
	flux_index=-2.0
	field=4.0
	l_min=0.0-0.5*field
	b_min=0.0-0.5*field
	num_sources=4
	
	# Loop over sources to be generated
	for i in range(0,num_sources):
	
		# Get random flux
		flux=math.log10(flux_min)+(math.log10(flux_max)-math.log10(flux_min))*random.random()
		flux=math.pow(10.0,flux)
	
		# Get random longitude/latitude
		l=l_min+field*random.random()
		b=b_min+field*random.random()
	
		# Add point-source Crabs
		sourcename="Source "+str(i)
		model.append(crab_like(sourcename, l, b, flux=flux, type="point"))
	
	# Return
	return model


# ================================================= #
# Setup a test point source for the survey analysis #
# ================================================= #
def add_test_sources(models, l=0.0, b=0.0, flux=1.0, index=-2.5):
	"""
	Add a test point-source with specified position, flux, and spectral index
	
	Arguments
	- l: longitude
	- b: latitude
	- flux: flux in Crab units
	- index: spectral index
	"""
	
	# Set source name
	sourcename="Test"
	
	# Add it to container
	models = add_source_model(models, any_pointlike(sourcename, l, b, \
	                          flux=flux, index=index, fitpos=False, fitidx=False))
	
	# Return
	return models


#==============================#
# Modified routine entry point #
#==============================#
if __name__ == '__main__':
	"""
	CTA survey simulation and analysis script.
	
	Notes:
	- Based on the initial script
	- Added functions: survey_square, add_survey_sources, add_random_sources
	- Added import: math, random
	
	"""
	# Set job parameters
	irf          = "kb_E_50h_v3"     # IRF
	survey_time  = 50.0*3600.0       # Total survey time
	pnt_step     = 1.0               # Spacing between pointings
	pnt_mintime  = 3600.0            # Minimum duration of a pointing
	pnt_fov      = 1.0               # Typical FoV size (useful to know where to probe sensitivity)
	pnt_maxnum   = 1                 # Maximum number of pointings
	emin         = 0.1               # minimum energy (TeV)
	emax         = 100.0             # maximum energy (TeV)
	#emin         = 1.0               # minimum energy (TeV)
	#emax         = 5.0               # maximum energy (TeV)
	nbins        = 20                # number of logarithmic energy bins
	test_thres   = 25.0              # Reference TS
	index        = -3.0              # Spectral index
	
	# Initialise response database
	caldb = "/usr/local/gamma/share/caldb/cta"	
	
	# Initialise flags
	need_help = False
	
	# Test for command line arguments
	if (len(sys.argv) > 1):
		surveyrange= float(sys.argv[1])
		
	if (len(sys.argv) > 2):
		surveystep= float(sys.argv[2])
		
	if (len(sys.argv) > 3):
		irf=sys.argv[3]			

	# Print help if needed and exit
	if need_help:
		print "Usage: example_survey.py [OPTIONS]"
		print "      galscan  Simulate Galactic plane scan"
		print "      extended Simulate also extended sources"
		print "     -h        Display this usage message"
		sys.exit()
	
	# Dump header
	print "***********************"
	print "* Simulate CTA survey *"
	print "***********************"
	
	# Remove any existing result files
	list = [glob.glob("*.fits"), glob.glob("*.log"), glob.glob("*.xml")]
	for files in list:
		for file in files:
			os.remove(file)
	
	# Set number of observations
	first        = True
	results_file = 'sensitivity.dat'
	for i in range(1,pnt_maxnum+1):
	
		# Initialise models
		bkg_model  = GModels()
		full_model = GModels()

		# Setup models
		bkg_model  = add_background_model(bkg_model)
		full_model = add_background_model(full_model)
		full_model = add_test_sources(full_model, 0.0, 0.0, flux=1.0, index=index)

		# Compute duration of pointings and starting longitude
		pnt_duration = survey_time/float(i)
		pnt_lmin     = 0.0-float(i-1)*pnt_step/2.0

		# Setup observations for survey 
		survey = survey_galplane(caldb, irf, lmin=pnt_lmin, lstep=pnt_step, \
							     nstep=i, duration=pnt_duration, \
							     emin=emin, emax=emax)
		
		print ""
		print "Simulating an array of "+str(i)+" observations along the plane"
	
		# Look for smallest flux for detection
		test_num   = 0
		test_max   = 50
		test_flux  = 0.001
		
		# Loop
		flux_value = []
		while True:
		
			# Update test counter
			test_num = test_num + 1
			
			# Set test flux
			full_model[1]['Prefactor'].value(test_flux)
			survey.models(full_model)
			
			# Simulate events
			obs = sim_obs(survey, nbins=nbins, seed=test_num)
		
			# Determine number of events in simulation
			nevents = 0.0
			for x in obs:
				nevents += x.events().number()

			# Fit background only
			obs.models(bkg_model)
			like = fit_obs(obs,log=False)
			
			# Store results
			LogL_bgm  = like.opt().value()
			Pre_bgm   = like.obs().models()[0]['Prefactor'].value()
			ePre_bgm  = like.obs().models()[0]['Prefactor'].error()
			Inx_bgm   = like.obs().models()[0]['Index'].value()
			eInx_bgm  = like.obs().models()[0]['Index'].error()
			npred_bgm = like.obs().npred()

			# Fit background and test source
			obs.models(full_model)
			like = fit_obs(obs,log=False)

			# Store results
			LogL_all  = like.opt().value()
			test_ts   = 2.0*(LogL_bgm-LogL_all)
			Pre_all   = like.obs().models()[0]['Prefactor'].value()
			ePre_all  = like.obs().models()[0]['Prefactor'].error()
			Inx_all   = like.obs().models()[0]['Index'].value()
			eInx_all  = like.obs().models()[0]['Index'].error()
			SPre_all  = like.obs().models()[1]['Prefactor'].value()
			eSPre_all = like.obs().models()[1]['Prefactor'].error()
			SInx_all  = like.obs().models()[1]['Index'].value()
			eSInx_all = like.obs().models()[1]['Index'].error()
			npred_all = like.obs().npred()
			if test_ts > 0:
				flux = math.sqrt(test_thres/test_ts)*SPre_all
			else:
				flux = 0.0
			flux_value.append(flux)
			
			# Compute sliding average
			flux = 0.0
			num  = 0.0
			for k in range(3):
				inx = test_num - k - 1
				if inx >= 0:
					flux += flux_value[inx]
					num  += 1.0
			flux /= num
			
			# Print results
			print "TS=%9.2f %8.5f lnL(BGD)=%8.1f lnL(ALL)=%8.1f #%5d" \
			      "  BGD: #%8.2f %5.2f+/-%4.2f %5.2f+/-%5.3f" \
			      "  ALL: #%8.2f %5.2f+/-%4.2f %5.2f+/-%5.3f" \
				  "  %8.5f+/-%7.5f %5.2f+/-%5.3f  Flux=%8.5f" % \
				  (test_ts, test_flux, LogL_bgm, LogL_all, nevents, \
				   npred_bgm, Pre_bgm, ePre_bgm, Inx_bgm, eInx_bgm, \
				   npred_all, Pre_all, ePre_all, Inx_all, eInx_all, \
				   SPre_all, eSPre_all, SInx_all, eSInx_all, \
				   flux)

			# Compare average flux to last average
			if test_num > 3:
				ratio = flux/test_flux
				if ratio >= 0.99 and ratio <= 1.01:
					print "Converged ("+str(ratio)+")"
					break
			
			# Use average for next iteration
			test_flux = flux
			
			# Exit loop if number of trials exhausted
			if (test_num >= test_max):
				print "Test ended after "+str(test_max)+" iterations."				
				break

		# Print results
		print "Sensitivity of "+str(i)+" observations: "+str(flux*1000.0)+" mCrab"
		
		# Write to CSV file
		result                 = {}
		colnames               = ['emin', 'emax', 'observations', 'sensitivity']
		result['emin']         = emin
		result['emax']         = emax
		result['observations'] = i
		result['sensitivity']  = flux
		if first:
			f      = open(results_file, 'w')
			writer = csv.DictWriter(f, colnames)
			writer.writerow(dict((_,_) for _ in colnames))
			first = False
		else:
			f = open(results_file, 'a')
		writer = csv.DictWriter(f, colnames)
		writer.writerow(result)
		f.close()
		
		# Delete some objects
		del obs
		del like

		
