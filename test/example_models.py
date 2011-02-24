#! /usr/bin/env python
# ===========================================================================================#
# This script illustrates how to build a unbinned analysis chain using CTAtools.
# ===========================================================================================#
from ctatools import *
from gammalib import *
from math import *
import os
import glob
import sys


# ================= #
# Analysis pipeline #
# ================= #
def pipeline(model_name):
	"""
	Unbinned analysis pipeline - keep intermediate results in memory.
    
    This function implements an analysis pipeline that successively calls
    ctobssim, ctselect and ctlike without saving the intermediate results as
	FITS files on disk. All data is only hold in memory.
	"""
	# Set script parameters
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

	# Write model name
	print "Model: "+model_name

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

	# Perform maximum likelihood fitting
	like = ctlike(select.obs())
	like.run()

	# Show model fitting results
	print like.obs().models()[0]
	
	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Perform unbinned analyses for various models.
	
	This script implements the analysis steps to perform an unbinned analysis.
	The steps include:
	- Event simulation
	- Event selection
	- Binned maximum likelihood fitting
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
		print "Usage: example_models.py [OPTIONS]"
		print "     -h       Display this usage message"
		print "     -noshow  Do not show data using matplotlib"
		sys.exit()

	# Dump header
	print "********************************"
	print "* CTA unbinned analysis script *"
	print "********************************"
		
	# Perform analysis for Crab model
	pipeline("data/crab.xml")

	# Perform analysis for disk model
	pipeline("data/disk.xml")

	# Perform analysis for Gaussian model
	pipeline("data/gauss.xml")
	
	# Perform analysis for shell model
	pipeline("data/shell.xml")
