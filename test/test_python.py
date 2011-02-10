#! /usr/bin/env python
# ===========================================================================================#
# This script illustrates how to build a unbinned analysis chain using CTAtools.
#
# If matplotlib is installed, the spectrum will be displayed on the screen. 
# ===========================================================================================#
from ctatools import *
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
	Unbinned analysis pipeline - save intermediate results.
    
    This function implements an analysis pipeline that successively calls
    ctobssim, ctselect and ctlike by saving the intermediate results as FITS
    files on disk. This replicates the way how the analysis would be done
    in a ftools-like approach.
	"""
	# Set script parameters
	model_name           = "data/crab.xml"
	events_name          = "events.fits"
	selected_events_name = "selected_events.fits"
	result_name          = "results.xml"
	caldb                = "irf"
	irf                  = "kb_E_50h_v3"
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
	sys.stdout.write(".")

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
	sys.stdout.write(".")

	# Perform maximum likelihood fitting
	like = ctlike()
	like.logFileOpen()
	like["evfile"].filename(selected_events_name)
	like["srcmdl"].filename(model_name)
	like["outmdl"].filename(result_name)
	like["method"].string("UNBINNED")
	like["caldb"].string(caldb)
	like["irf"].string(irf)
	like.execute()
	sys.stdout.write(".")

	# Return
	return


# ============================= #
# Analysis pipeline - version 2 #
# ============================= #
def pipeline_v2():
	"""
	Unbinned analysis pipeline - keep intermediate results in memory.
    
    This function implements an analysis pipeline that successively calls
    ctobssim, ctselect and ctlike without saving the intermediate results as
	FITS files on disk. All data is only hold in memory.
	"""
	# Set script parameters
	model_name  = "data/crab.xml"
	caldb       = "irf"
	irf         = "kb_E_50h_v3"
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
	sys.stdout.write(".")

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
	sys.stdout.write(".")

	# Perform maximum likelihood fitting
	like = ctlike(select.obs())
	like.run()
	sys.stdout.write(".")

	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Test ctatools Python modules.
	"""
	# Dump header
	print ""
	print "*************************************"
	print "* ctatools Python interface testing *"
	print "************************************"
	
	# Remove any existing result files
	list = [glob.glob("*.fits"), glob.glob("*.log"), glob.glob("*.xml")]
	for files in list:
		for file in files:
			os.remove(file)
	
	# Save intermediate results on disk
	print "Test executable analysis: ",
	pipeline_v1()
	print " ok"

	# Keep results in memory
	print "Test in-memory analysis: ",
	pipeline_v2()
	print " ok"
