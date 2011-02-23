#! /usr/bin/env python
# ===========================================================================================#
# This script illustrates how to build a binned analysis pipeline using CTAtools.
#
# If matplotlib is installed, a counts spectrum and an offset histogram will be displayed. 
# ===========================================================================================#
from ctatools import *
from gammalib import *
from math import *
import os
import glob
import sys
import csv

# Parameters
events_name = "events.fits"
cntmap_name = "cntmap.fits"
result_name = "results.xml"
caldb       = "irf"
irf         = "kb_E_50h_v3"
ra          =   83.63
dec         =   22.01
rad_sim     =   10.0
tstart      =    0.0
tstop       =    1e4
emin        =    1
emax        =   10
enumbins    =    5
nxpix       =  200
nypix       =  200
binsz       =    0.02
coordsys    = "CEL"
proj        = "CAR"

def run_sim(seed, model_name):
	sim = ctobssim()
	sim.logFileOpen()  # We need this to explicitely open the log file in Python mode
	sim["infile"].filename(model_name)
	sim["outfile"].filename(events_name)
	sim["caldb"].string(caldb)
	sim["irf"].string(irf)
	sim["seed"].integer(seed)
	sim["ra"].real(ra)
	sim["dec"].real(dec)
	sim["rad"].real(rad_sim)
	sim["tmin"].real(tstart)
	sim["tmax"].real(tstop)
	sim["emin"].real(emin)
	sim["emax"].real(emax)
	sim.execute()
	return sim

def run_bin():
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
	return bin
	
def run_like(model_name):
	like = ctlike()
	like.logFileOpen()  # We need this to explicitely open the log file in Python mode
	like["cntmap"].filename(cntmap_name)
	like["srcmdl"].filename(model_name)
	like["outmdl"].filename(result_name)
	like["method"].string("BINNED")
	like["caldb"].string(caldb)
	like["irf"].string(irf)
	like.execute()
	return like

def test_coverage(n_trials=3, model_name="data/crab.xml"):
	"""Simulate and fit ntrials time and look at the
	pull distribution for all free parameters.
	It should be Gaussian with mean 0 and rms 1"""
	
	results = []
	colnames = [] # to preserve the order!
	for seed in range(n_trials):
		print 'seed:', seed
		sim = run_sim(seed, model_name)
		run_bin()
		like = run_like(model_name)
		models = like.obs().models()
		result = {}
		colnames = []
		for par in range(models.npars()):
			p = models.par(par)
			if p.isfree():
				name = p.name()
				colnames.append(name)
				fitted_value = p.value()
				real_value = sim.obs().models().par(par).value()
				error = p.error()
				pull = (fitted_value - real_value) / error
				result[name] = pull
		results.append(result)
	return results, colnames

if __name__ == '__main__':
	usage = "test_coverage xml_model_file n_trials results_file"
	if len(sys.argv) < 4:
		print usage
		sys.exit()
	
	xml_model_file = sys.argv[1]
	n_trials       = int(sys.argv[2])
	results_file   = sys.argv[3]
	
	print("Testing coverage for file {0} with {1} trials"
		"".format(xml_model_file, n_trials))
	results, colnames = test_coverage(n_trials, xml_model_file)
	
	print("Writing pull values for all free parameters to file {0}"
		"".format(results_file))
	
	f = open(results_file, 'w')
	writer = csv.DictWriter(f, colnames)
	writer.writerow(dict((_,_) for _ in colnames))
	writer.writerows(results)
	f.close()