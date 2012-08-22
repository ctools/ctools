#! /usr/bin/env python
# ==========================================================================
# This script illustrates how to build a binned analysis pipeline using the
# ctools. If matplotlib is installed, a counts spectrum and an offset
# histogram will be displayed.
#
# Required 3rd party modules:
# None
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
from math import *
import os
import glob
import sys


# ============================= #
# Analysis pipeline - version 1 #
# ============================= #
def pipeline_v1():
	"""
	Binned analysis pipeline - save intermediate results.
    
    This function implements an analysis pipeline that successively calls
    ctobssim, ctbin and ctlike by saving the intermediate results as FITS
    files on disk. This replicates the way how the analysis would be done
    in a ftools-like approach.
	"""
	# Set script parameters
	model_name  = "${GAMMALIB}/share/models/crab.xml"
	events_name = "events.fits"
	cntmap_name = "cntmap.fits"
	result_name = "results.xml"
	caldb       = "${GAMMALIB}/share/caldb/cta"
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
	sim.logFileOpen()  # We need this to explicitely open the log file in Python mode
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
	sim.execute()
	sys.stdout.write("Simulated events ("+str(sim.celapse())+" CPU seconds)\n")

	# Update timing
	wall_seconds += sim.telapse()
	cpu_seconds  += sim.celapse()
	
	# Bin events into counts map
	bin = ctbin()
	bin.logFileOpen()  # We need this to explicitely open the log file in Python mode
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
	bin.execute()
	sys.stdout.write("Binned events into counts map ("+str(bin.celapse())+" CPU seconds)\n")

	# Update timing
	wall_seconds += bin.telapse()
	cpu_seconds  += bin.celapse()
	
	# Perform maximum likelihood fitting
	like = ctlike()
	like.logFileOpen()  # We need this to explicitely open the log file in Python mode
	like["infile"].filename(cntmap_name)
	like["srcmdl"].filename(model_name)
	like["outmdl"].filename(result_name)
	like["caldb"].string(caldb)
	like["irf"].string(irf)
	like.execute()
	sys.stdout.write("Maximum likelihood fitting ("+str(like.celapse())+" CPU seconds)\n")

	# Update timing
	wall_seconds += like.telapse()
	cpu_seconds  += like.celapse()
	
	# Show total times
	sys.stdout.write("Total wall time elapsed: "+str(wall_seconds)+" seconds\n")
	sys.stdout.write("Total CPU time used ...: "+str(cpu_seconds)+" seconds\n")

	# Return
	return


# ============================= #
# Analysis pipeline - version 2 #
# ============================= #
def pipeline_v2():
	"""
	Binned analysis pipeline - keep intermediate results in memory.
    
    This function implements an analysis pipeline that successively calls
    ctobssim, ctbin and ctlike without saving the intermediate results as
	FITS files on disk. All data are hold in memory.
	
	At the end, results are plotted (if matplotlib is installed)
	"""
	# Set script parameters
	model_name  = "${GAMMALIB}/share/models/crab.xml"
	caldb       = "${GAMMALIB}/share/caldb/cta"
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
	sys.stdout.write("Simulated events ("+str(sim.celapse())+" CPU seconds)\n")

	# Update timing
	wall_seconds += sim.telapse()
	cpu_seconds  += sim.celapse()

	# Bin events into counts map
	bin = ctbin(sim.obs())
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
	sys.stdout.write("Binned events into counts map ("+str(bin.celapse())+" CPU seconds)\n")

	# Update timing
	wall_seconds += bin.telapse()
	cpu_seconds  += bin.celapse()

	# Perform maximum likelihood fitting
	like = ctlike(bin.obs())
	like.run()
	sys.stdout.write("Maximum likelihood fitting ("+str(like.celapse())+" CPU seconds)\n")

	# Update timing
	wall_seconds += like.telapse()
	cpu_seconds  += like.celapse()
		
	# Show total times
	sys.stdout.write("Total wall time elapsed: "+str(wall_seconds)+" seconds\n")
	sys.stdout.write("Total CPU time used ...: "+str(cpu_seconds)+" seconds\n")

	# Show model fitting results
	#sys.stdout.write(like.obs().models()+"\n")

	# Plot counts
	plot_counts(bin.obs())
	
	# Return
	return


# =========== #
# Plot counts #
# =========== #
def plot_counts(observations):
	"""
	Plot counts.
	"""
	# Only proceed if matplotlib is available
	try:
		# Import matplotlib
		import matplotlib.pyplot as plt

		# Set legend fontsize
		params = {'legend.fontsize': 10}
		plt.rcParams.update(params)

		# Set plot styles
		styles = ['b-', 'g-', 'y-', 'n-']

		# Dump header
		sys.stdout.write("\n")
		sys.stdout.write("Make plots (using matplotlib):\n")
		sys.stdout.write("==============================\n")
		
		# Create figure 1
		plt.figure(1,figsize=(12,6))
		plt.subplots_adjust(hspace=.7)

		# Create subplot 1
		plt.subplot(121)
		plt.title("Spectrum (summed over all pixels)")
		
		# Loop over observations
		for obs in observations:
		
			# Get event cube
			cube = cast_GCTAEventCube(obs.events())
			
			# Create energy axis
			energy  = []
			ebounds = cube.ebounds()
			for i in range(ebounds.size()):
				energy.append(ebounds.elogmean(i).TeV())

			# Create spectrum
			sys.stdout.write("Extract data:\n")
			counts = [0.0 for i in range(ebounds.size())]
			for bin in cube:
				index         = ebounds.index(bin.energy())
				counts[index] = counts[index] + bin.counts()

			# Create error bars
			error = [sqrt(c) for c in counts]

			# Plot spectrum
			plt.loglog(energy, counts, 'ro', label='data')
			plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

			# Extract models
			sys.stdout.write("Extract models:\n")
			sum_model = [0.0 for i in range(ebounds.size())]
			for k, m in enumerate(observations.models()):
				sys.stdout.write("- "+m.name()+"\n")
				model = [0.0 for i in range(ebounds.size())]
				for bin in cube:
					index        = ebounds.index(bin.energy())
					prob         = m.eval(bin, obs)
					model[index] = model[index] + prob * bin.size()
				for i in range(ebounds.size()):
					sum_model[i] = sum_model[i] + model[i]
				plt.loglog(energy, model, styles[k], label=m.name())
			plt.loglog(energy, sum_model, 'r-', label='total')

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
		for obs in observations:
		
			# Get event cube
			cube = cast_GCTAEventCube(obs.events())
			
			# Create offset histogram
			sys.stdout.write("Extract data:\n")
			nx       = 30
			doffset2 = 0.01
			offset2  = [(i+0.5)*doffset2 for i in range(nx)]
			counts   = [0.0  for i in range(nx)]
			for bin in cube:
				off   = bin.dir().dist_deg(crab)
				off2  = off*off
				index = int(off2/doffset2)
				if index < nx:
					counts[index] = counts[index] + bin.counts()

			# Create error bars
			error = [sqrt(c) for c in counts]

			# Plot distribution
			plt.semilogy(offset2, counts, 'ro', label='data')

			# Extract models
			sys.stdout.write("Extract models:\n")
			sum_model = [0.0 for i in range(nx)]
			for k, m in enumerate(observations.models()):
				sys.stdout.write("- "+m.name()+"\n")
				model = [0.0 for i in range(nx)]
				for bin in cube:
					off   = bin.dir().dist_deg(crab)
					off2  = off*off
					index = int(off2/doffset2)
					if index < nx:
						prob         = m.eval(bin, obs)
						model[index] = model[index] + prob * bin.size()
				for i in range(nx):
					sum_model[i] = sum_model[i] + model[i]
				plt.plot(offset2, model, styles[k], label=m.name())
			plt.plot(offset2, sum_model, 'r-', label='total')
			plt.ylim(ymin=0.1)

		# Put labels
		plt.xlabel("Offset (deg^2)")
		plt.ylabel("Counts")
		plt.legend(loc="upper right")

		# Show counts spectra
		plt.show()

	except:
		sys.stdout.write("Matplotlib is not (correctly) installed on your system. No counts spectra are shown.\n")

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
	
	Two variations of the script are implemented. One that calls the classes
	like executables, storing the intermediate results on disk, and a second
	that uses the results of the precedent class directly without storing
	any data on disk.
	
	The second variation will also show a counts spectrum of the data after
	fitting (only if matplotlib is available on your system). Note that the
	counts spectrum is summed over the entire data space, i.e. no event selection
	around the Crab has been performed.
	"""
	# Dump header
	sys.stdout.write("*************************************\n")
	sys.stdout.write("*    CTA binned analysis scripts    *\n")
	sys.stdout.write("*************************************\n")
	sys.stdout.write("... this script will take 1-2 minutes\n")

	# Remove any existing result files
	list = [glob.glob("*.fits"), glob.glob("*.log"), glob.glob("*.xml")]
	for files in list:
		for file in files:
			os.remove(file)
	
	# Save intermediate results on disk
	sys.stdout.write("\n")
	sys.stdout.write("Executable analysis pipeline:\n")
	sys.stdout.write("=============================\n")
	pipeline_v1()

	# Do analysis in memory
	sys.stdout.write("\n")
	sys.stdout.write("In memory analysis pipeline:\n")
	sys.stdout.write("============================\n")
	pipeline_v2()
	sys.stdout.write("\n")
