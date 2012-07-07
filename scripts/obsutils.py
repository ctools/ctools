# ==========================================================================
# This script provides a number of functions that are useful for handling
# CTA observations.
#
# Copyright (C) 2011-2012 Juergen Knoedlseder
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


# ===================== #
# Simulate observations #
# ===================== #
def sim(obs, log=False, debug=False, seed=0, nbins=0, binsz=0.05, npix=200):
	"""
	Simulate events for all observations in the container.
	
	Parameters:
	 obs   - Observation container
	Keywords:
	 log   - Create log file(s)
	 debug - Create screen dump
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

	# Optionally switch-on debugging model
	if debug:
		sim["debug"].boolean(True)
	
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

		# Optionally switch-on debugging model
		if debug:
			bin["debug"].boolean(True)

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
    
	# Delete the simulation
	del sim
	
	# Return observation container
	return obs


# ================ #
# Fit observations #
# ================ #
def fit(obs, log=False, debug=False):
	"""
	Perform maximum likelihood fitting of observations in the container.
	
	Parameters:
	 obs   - Observation container
	Keywords:
	 log   - Create log file(s)
	 debug - Create screen dump
	"""
	# Allocate ctlike application
	like = ctlike(obs)
	
	# Optionally open the log file
	if log:
		like.logFileOpen()
	
	# Optionally switch-on debugging model
	if debug:
		like["debug"].boolean(True)

	# Run ctlike application.
	like.run()
	
	# Return observations
	return like


# ================= #
# Create counts map #
# ================= #
def cntmap(obs, proj="TAN", coord="GAL", xval=0.0, yval=0.0, \
           binsz=0.05, nxpix=200, nypix=200, \
		   outname="cntmap.fits"):
	"""
	Creates a counts map by combining the events of all observations.
	The counts map will be a summed map over all energies.
	
	Parameters:
	 obs     - Observation container
	Keywords:
	 proj    - Projection type (e.g. TAN, CAR, STG, ...) (default: TAN)
	 coord   - Coordinate type (GAL, CEL) (default: GAL)
	 xval    - Reference longitude value [deg] (default: 0.0)
	 yval    - Reference latitude value [deg] (default: 0.0)
	 binsz   - Pixel size [deg/pixel] (default: 0.05)
	 nxpix   - Number of pixels in X direction (default: 200)
	 nypix   - Number of pixels in Y direction (default: 200)
	 outname - Counts map FITS filename (default: cntmap.fits)
	"""
	# Allocate counts map
	map = GSkymap(proj, coord, xval, yval, -binsz, binsz, nxpix, nypix, 1)
	
	# Fill all observations
	for run in obs:
		
		# Loop over all events
		for event in run.events():
			
			# Cast to CTA event
			event = cast_GCTAEventAtom(event)
			
			# Determine sky pixel
			skydir = event.dir().skydir()
			pixel  = map.dir2pix(skydir)
			
			# Set pixel
			map[pixel] += 1.0
	
	# Save sky map. The clobber flag is set to True, so any existing FITS
	# file will be overwritten.
	map.save(outname, True)
	
	# Return counts map
	return map


# ================ #
# Create model map #
# ================ #
def modmap(obs, eref=0.1, proj="TAN", coord="GAL", xval=0.0, yval=0.0, \
           binsz=0.05, nxpix=200, nypix=200, \
           outname="modmap.fits"):
	"""
	Make model map for a given reference energy by combining all observations.
	The model map will be evaluated for a given reference energy 'eref' and will
	be given in units of [counts/(sr MeV s)].
	
	Parameters:
	 obs     - Observation container
	Keywords:
	 eref    - Reference energy for which model is created [TeV] (default: 0.1)
	 proj    - Projection type (e.g. TAN, CAR, STG, ...) (default: TAN)
	 coord   - Coordinate type (GAL, CEL) (default: GAL)
	 xval    - Reference longitude value [deg] (default: 0.0)
	 yval    - Reference latitude value [deg] (default: 0.0)
	 binsz   - Pixel size [deg/pixel] (default: 0.05)
	 nxpix   - Number of pixels in X direction (default: 200)
	 nypix   - Number of pixels in Y direction (default: 200)
	 outname - Model map FITS filename (default: modmap.fits)
	"""
	# Allocate model map
	map = GSkymap(proj, coord, xval, yval, -binsz, binsz, nxpix, nypix, 1)
	
	# Set reference energy, time and direction. The time is not initialised and is
	# in fact not used (as the IRF is assumed to be time independent for now).
	# The sky direction is set later using the pixel values.
	energy  = GEnergy()
	time    = GTime()
	instdir = GCTAInstDir()
	energy.TeV(eref)
	
	# Loop over all map pixels
	for pixel in range(map.npix()):
		
		# Get sky direction
		skydir = map.pix2dir(pixel)
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
		map[pixel] = value
	
	# Save sky map
	map.save(outname, True)
	
	# Return model map
	return map


# ======================= #
# Set one CTA observation #
# ======================= #
def set(pntdir, tstart=0.0, duration=1800.0, deadc=0.95, \
        emin=0.1, emax=100.0, rad=5.0, \
        irf="kb_E_50h_v3", caldb="$GAMMALIB/share/caldb/cta"):
	"""
	Returns a single CTA observation. By looping over this function we can
	add CTA observations to the observation container.
	
	Parameters:
	 pntdir   - Pointing direction
	Keywords:
	 tstart   - Start time [seconds] (default: 0.0)
	 duration - Duration of observation [seconds] (default: 1800.0)
	 deadc    - Deadtime correction factor (default: 0.95)
	 emin     - Minimum event energy [TeV] (default: 0.1)
	 emax     - Maximum event energy [TeV] (default: 100.0)
	 rad      - ROI radius used for analysis [deg] (default: 5.0)
	 irf      - Instrument response function (default: kb_E_50h_v3)
	 caldb    - Calibration database path (default: $GAMMALIB/share/caldb/cta)
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
	
	# Set ontime, livetime, and deadtime correction factor
	obs.ontime(duration)
	obs.livetime(duration*deadc)
	obs.deadc(deadc)
	
	# Return observation
	return obs
