#! /usr/bin/env python
# ===========================================================================================#
# This script illustrates how to build a binned analysis chain using CTAtools.
#
# If matplotlib is installed, the spectrum will be displayed on the screen. 
# ===========================================================================================#
from ctatools import *
from math import *


# =================== #
# Processing pipeline #
# =================== #
def process():
	"""
	Binned analysis processing pipeline.
	"""
	# Set script parameters
	model_name  = "data/crab.xml"
	events_name = "events.fits"
	cntmap_name = "cntmap.fits"
	result_name = "results.xml"
	caldb       = "irf"
	irf         = "kb_E_50h_v3"
	ra          =   83.63
	dec         =   22.01
	rad_sim     =   10.0
	tstart      =    0.0
	tstop       = 1800.0
	emin        =    0.1
	emax        =  100.0
	enumbins    =   20
	nxpix       =  200
	nypix       =  200
	binsz       =    0.02
	coordsys    = "CEL"
	proj        = "CAR"

	# Initialise timing
	wall_seconds = 0.0
	cpu_seconds  = 0.0

	# Simulate events
	sim = ctobssim()
	sim["infile"].filename(model_name)
	sim["outfile"].filename(events_name)
	sim["caldb"].string(caldb)
	sim["irf"].string(irf)
	sim["ra"].real(ra)
	sim["dec"].real(dec)
	sim["rad"].real(rad_sim)
	sim["tmin"].real(tstart)
	sim["tmax"].real(tstop)
	sim["emin"].real(emin)
	sim["emax"].real(emax)
	sim.run()
	print "Simulated events ("+str(sim.celapse())+" CPU seconds)"

	# Update timing
	wall_seconds += sim.telapse()
	cpu_seconds  += sim.celapse()
	
	# Bin events into counts map
	bin = ctbin()
	bin["evfile"].filename(events_name)
	bin["outfile"].filename(cntmap_name)
	bin["emin"].real(emin)
	bin["emax"].real(emax)
	bin["enumbins"].integer(enumbins)
	bin["nxpix"].integer(nxpix)
	bin["nypix"].integer(nypix)
	bin["binsz"].real(binsz)
	bin["coordsys"].string(coordsys)
	bin["xref"].real(ra)
	bin["yref"].real(dec)
	bin["proj"].string(proj)
	bin.run()
	print "Binned events into counts map ("+str(bin.celapse())+" CPU seconds)"

	# Update timing
	wall_seconds += bin.telapse()
	cpu_seconds  += bin.celapse()
	
	# Perform maximum likelihood fitting
	like = ctlike()
	like["cntmap"].filename(cntmap_name)
	like["srcmdl"].filename(model_name)
	like["outmdl"].filename(result_name)
	like["method"].string("BINNED")
	like["stat"].string("POISSON")
	like["caldb"].string(caldb)
	like["irf"].string(irf)
	like.run()
	print "Maximum likelihood fitting ("+str(like.celapse())+" CPU seconds)"

	# Update timing
	wall_seconds += like.telapse()
	cpu_seconds  += like.celapse()
	
	# Show total times
	print "Total wall time elapsed: "+str(wall_seconds)+" seconds"
	print "Total CPU time used ...: "+str(cpu_seconds)+" seconds"

	# Return
	return


# ============ #
# Show photons #
# ============ #
def show_photons(photons, xmlname, e_min, e_max, area, duration, ebins=30):
	"""
	Show photons using matplotlib (if available).
	"""
	# Only proceed if matplotlib is available
	try:
		# Import matplotlib
		import matplotlib.pyplot as plt

		# Create figure
		plt.figure(1)
		plt.title("MC simulated photon spectrum ("+str(e_min)+'-'+str(e_max)+" TeV)")
		
		# Setup energy range covered by data
		emin = GEnergy()
		emax = GEnergy()
		emin.TeV(e_min)
		emax.TeV(e_max)
		ebds = GEbounds()
		ebds.setlog(emin, emax, ebins)
		
		# Create energy axis
		energy = []
		for i in range(ebds.size()):
			energy.append(ebds.elogmean(i).TeV())

		# Fill histogram
		counts = [0.0 for i in range(ebds.size())]
		for photon in photons:
			index         = ebds.index(photon.energy())
			counts[index] = counts[index] + 1.0

		# Create error bars
		error = [sqrt(c) for c in counts]

		# Get model values
		models = GModels(xmlname)
		crab   = cast_GModelSky(models[0])
		model = []
		d = GSkyDir()
		d.radec_deg(83.6331, 22.0145)
		t = GTime()
		for i in range(ebds.size()):
			eval   = ebds.elogmean(i)
			ewidth = ebds.emax(i) - ebds.emin(i)
			f      = crab.value(d, eval, t) * area * duration * ewidth.MeV()
			model.append(f)

		# Plot data
		plt.loglog(energy, counts, 'ro')
		plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

		# Plot model
		plt.plot(energy, model, 'b-')

		# Set axes
		plt.xlabel("Energy (TeV)")
		plt.ylabel("Number of incident photons")

		# Allocate histogram
		# Show plot
		plt.show()

	except:
		print "Matplotlib is not (correctly) installed on your system. No data are shown."

	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	CTA binned analysis script.
	
	This script implements the analysis steps to perform a binned analysis.
	The steps include:
	- Event simulation
	- Event binning
	- Binned maximum likelihood fitting
	"""
	# Dump header
	print "*****************************"
	print "* CTA binned analysis script *"
	print "*****************************"

	# Process data
	process()
	
	# Show photons
	#show_photons(photons, xmlname, e_min, e_max, area, duration)
