#! /usr/bin/env python
# ==========================================================================
# This script illustrates how to build an unbinned analysis pipeline with
# the ctools. If matplotlib is installed, the spectrum will be displayed
# on the screen.
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
try:
	import matplotlib.pyplot as plt
	has_matplotlib = True
except:
	has_matplotlib = False


# ============================= #
# Analysis pipeline - version 1 #
# ============================= #
def pipeline_v1():
	"""
	Unbinned analysis pipeline - save intermediate results.
    
    This function implements an analysis pipeline that successively calls
    ctobssim, ctselect and ctlike by saving the intermediate results as FITS
    files on disk. This replicates the way how the analysis would be done
    in a ftools-like approach.
	"""
	# Set script parameters
	model_name           = "${CTOOLS}/share/models/crab.xml"
	events_name          = "events.fits"
	selected_events_name = "selected_events.fits"
	result_name          = "results.xml"
	caldb                = "${CALDB}/data/cta/dummy/bcf"
	irf                  = "cta_dummy_irf"
	ra                   =   83.63
	dec                  =   22.01
	rad_sim              =   10.0
	tstart               =    0.0
	tstop                = 1800.0
	emin                 =    0.1
	emax                 =  100.0
	rad_select           =    3.0

	# Initialise timing
	wall_seconds = 0.0
	cpu_seconds  = 0.0

	# Simulate events
	sim = ctobssim()
	sim.logFileOpen()
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
	print("Simulated events ("+str(sim.celapse())+" CPU seconds)")
	
	# Update timing
	wall_seconds += sim.telapse()
	cpu_seconds  += sim.celapse()
	
	# Select events
	select = ctselect()
	select.logFileOpen()
	select["infile"].filename(events_name)
	select["outfile"].filename(selected_events_name)
	select["ra"].real(ra)
	select["dec"].real(dec)
	select["rad"].real(rad_select)
	select["tmin"].real(tstart)
	select["tmax"].real(tstop)
	select["emin"].real(emin)
	select["emax"].real(emax)
	select.execute()
	print("Select event ("+str(select.celapse())+" CPU seconds)")

	# Update timing
	wall_seconds += select.telapse()
	cpu_seconds  += select.celapse()
	
	# Perform maximum likelihood fitting
	like = ctlike()
	like.logFileOpen()
	like["infile"].filename(selected_events_name)
	like["srcmdl"].filename(model_name)
	like["outmdl"].filename(result_name)
	like["caldb"].string(caldb)
	like["irf"].string(irf)
	like.execute()
	print("Maximum likelihood fitting ("+str(like.celapse())+" CPU seconds)")

	# Update timing
	wall_seconds += like.telapse()
	cpu_seconds  += like.celapse()
	
	# Show total times
	print("Total wall time elapsed: "+str(wall_seconds)+" seconds")
	print("Total CPU time used ...: "+str(cpu_seconds)+" seconds")

	# Return
	return


# ============================= #
# Analysis pipeline - version 2 #
# ============================= #
def pipeline_v2(show_data):
	"""
	Unbinned analysis pipeline - keep intermediate results in memory.
    
    This function implements an analysis pipeline that successively calls
    ctobssim, ctselect and ctlike without saving the intermediate results as
	FITS files on disk. All data is only hold in memory.
	"""
	# Set script parameters
	model_name  = "${CTOOLS}/share/models/crab.xml"
	caldb       = "${CALDB}/data/cta/dummy/bcf"
	irf         = "cta_dummy_irf"
	ra          =   83.63
	dec         =   22.01
	rad_sim     =   10.0
	tstart      =    0.0
	tstop       = 1800.0
	emin        =    0.1
	emax        =  100.0
	rad_select  =    3.0

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
	print("Simulated events ("+str(sim.celapse())+" CPU seconds)")

	# Update timing
	wall_seconds += sim.telapse()
	cpu_seconds  += sim.celapse()

	# Select events
	select = ctselect(sim.obs())
	select["ra"].real(ra)
	select["dec"].real(dec)
	select["rad"].real(rad_select)
	select["tmin"].real(tstart)
	select["tmax"].real(tstop)
	select["emin"].real(emin)
	select["emax"].real(emax)
	select.run()
	print("Select event ("+str(select.celapse())+" CPU seconds)")

	# Update timing
	wall_seconds += select.telapse()
	cpu_seconds  += select.celapse()

	# Perform maximum likelihood fitting
	like = ctlike(select.obs())
	like.run()
	print("Maximum likelihood fitting ("+str(like.celapse())+" CPU seconds)")

	# Update timing
	wall_seconds += like.telapse()
	cpu_seconds  += like.celapse()
		
	# Show total times
	print("Total wall time elapsed: "+str(wall_seconds)+" seconds")
	print("Total CPU time used ...: "+str(cpu_seconds)+" seconds")

	# Show model fitting results
	#print like.obs().models()
	
	# Optionally plot counts
	if show_data:
		if has_matplotlib:
			plot_counts(like.obs())
		else:
			sys.stdout.write("Matplotlib is not (correctly) installed on your system. No counts spectra are shown.\n")
	
	# Return
	return


# =========== #
# Plot counts #
# =========== #
def plot_counts(observations):
	"""
	Plot counts.
	"""
	# Set legend fontsize
	params = {'legend.fontsize': 10}
	plt.rcParams.update(params)

	# Set plot styles
	styles = ['b-', 'g-', 'y-', 'n-']

	# Dump header
	print("")
	print("Make plots (using matplotlib):")
	print("==============================")
	
	# Create figure 1
	plt.figure(1,figsize=(12,6))
	plt.subplots_adjust(hspace=.7)

	# Create subplot 1
	plt.subplot(121)
	plt.title("Spectrum (summed over all pixels)")
	
	# Loop over observations
	for obs in observations:
	
		# Get event list
		list = obs.events()
		
		# Create energy axis
		ebounds = GEbounds(10, GEnergy(0.1, "TeV"), GEnergy(100.0, "TeV"))
		energy  = [ebounds.elogmean(i).TeV() for i in range(ebounds.size())]

		# Create spectrum
		print("Extract data:")
		counts = [0.0 for i in range(ebounds.size())]
		print(list.size())
		for atom in list:
			index         = ebounds.index(atom.energy())
			counts[index] = counts[index] + 1.0

		# Create error bars
		error = [sqrt(c) for c in counts]

		# Plot spectrum
		plt.loglog(energy, counts, 'ro', label='data')
		#plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

		# Extract models
		#print("Extract models:")
		#sum_model = [0.0 for i in range(ebounds.size())]
		#for k, m in enumerate(observations.models()):
		#	print("- "+m.name())
		#	model = [0.0 for i in range(ebounds.size())]
		#	for atom in list:
		#		index        = ebounds.index(atom.energy())
		#		prob         = m.eval(atom, obs)
		#		model[index] = model[index] + prob * atom.size()
		#	for i in range(ebounds.size()):
		#		sum_model[i] = sum_model[i] + model[i]
		#	plt.loglog(energy, model, styles[k], label=m.name())
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
	for obs in observations:
	
		# Get event cube
		list = obs.events()
		
		# Create offset histogram
		print "Extract data:"
		nx       = 30
		doffset2 = 0.01
		offset2  = [(i+0.5)*doffset2 for i in range(nx)]
		counts   = [0.0  for i in range(nx)]
		for atom in list:
			off   = GCTAInstDir(atom.dir()).dir().dist_deg(crab)
			off2  = off*off
			index = int(off2/doffset2)
			if index < nx:
				counts[index] = counts[index] + 1.0

		# Create error bars
		error = [sqrt(c) for c in counts]

		# Plot distribution
		plt.semilogy(offset2, counts, 'ro', label='data')

		# Extract models
		print "Extract models:"
		sum_model = [0.0 for i in range(nx)]
		for k, m in enumerate(observations.models()):
			print "- "+m.name()
			model = [0.0 for i in range(nx)]
			for atom in list:
				off   = GCTAInstDir(atom.dir()).dir().dist_deg(crab)
				off2  = off*off
				index = int(off2/doffset2)
				if index < nx:
					prob         = m.eval(atom, obs)
					model[index] = model[index] + prob * atom.size()
			for i in range(nx):
				sum_model[i] = sum_model[i] + model[i]
			#plt.plot(offset2, model, styles[k], label=m.name())
		#plt.plot(offset2, sum_model, 'r-', label='total')
		#plt.ylim(ymin=0.1)

	# Put labels
	plt.xlabel("Offset (deg^2)")
	plt.ylabel("Counts")
	plt.legend(loc="upper right")

	# Show counts spectra
	plt.show()

	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	CTA unbinned analysis script.
	
	This script implements the analysis steps to perform an unbinned analysis.
	The steps include:
	- Event simulation
	- Event selection
	- Binned maximum likelihood fitting
	
	Two variations of the script are implemented. One that calls the classes
	like executables, storing the intermediate results on disk, and a second
	that uses the results of the precedent class directly without storing
	any data on disk.
	"""
	# Initialise flags
	need_help = False
	show_data = True
	
	# Test for command line arguments
	print sys.argv[0]
	if (len(sys.argv) > 1):
		if sys.argv[1] == "-h":
			need_help = True
		elif sys.argv[1] == "-noshow":
			show_data = False
		else:
			need_help = True			

	# Print help if needed and exit
	if need_help:
		print "Usage: example_binned.py [OPTIONS]"
		print "     -h       Display this usage message"
		print "     -noshow  Do not show data using matplotlib"
		sys.exit()

	# Dump header
	print "*********************************"
	print "* CTA unbinned analysis scripts *"
	print "*********************************"
	
	# Remove any existing result files
	list = [glob.glob("*.fits"), glob.glob("*.log"), glob.glob("*.xml")]
	for files in list:
		for file in files:
			os.remove(file)
	
	# Save intermediate results on disk
	print ""
	print "Executable analysis pipeline:"
	print "============================="
	pipeline_v1()
	
	# Do analysis in memory
	print ""
	print "In memory analysis pipeline:"
	print "============================"
	pipeline_v2(show_data)
	print ""
