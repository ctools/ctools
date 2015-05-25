#! /usr/bin/env python
# ==========================================================================
# This script performs an unbinned maximum likelihood analysis for a
# variety of models. This allows checking the models.
#
# Usage:
#   ./check_models.py
#
# ==========================================================================
#
# Copyright (C) 2011-2015 Juergen Knoedlseder
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
import gammalib
import ctools
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
	caldb       = "prod2"
	irf         = "South_50h"
	ra          =   83.63
	dec         =   22.01
	rad_sim     =   10.0
	tstart      =    0.0
	tstop       = 1800.0
	emin        =    0.1
	emax        =  100.0
	rad_select  =    3.0

	# Write model name
	print("*** Model: "+model_name+" ************************************")

	# Initialise timing
	wall_seconds = 0.0
	cpu_seconds  = 0.0

	# Simulate events
	sim = ctools.ctobssim()
	sim["inmodel"].filename(model_name)
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
	select = ctools.ctselect(sim.obs())
	select["ra"].real(ra)
	select["dec"].real(dec)
	select["rad"].real(rad_select)
	select["tmin"].real(tstart)
	select["tmax"].real(tstop)
	select["emin"].real(emin)
	select["emax"].real(emax)
	select.run()

	# Perform maximum likelihood fitting
	like = ctools.ctlike(select.obs())
	like.run()

	# Show model fitting results
	print(like.obs().models()[0])
	
	# Return
	return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
	"""
	Perform unbinned analyses for various models.
	"""
	# Initialise flags
	need_help = False
	
	# Test for command line arguments
	print(sys.argv[0])
	if (len(sys.argv) > 1):
		if sys.argv[1] == "-h":
			need_help = True
		else:
			need_help = True			

	# Print help if needed and exit
	if need_help:
		print("Usage: check_models.py [OPTIONS]")
		print("     -h       Display this usage message")
		sys.exit()

	# Dump header
	print("*******************************************")
	print("* Check models using an unbinned analysis *")
	print("*******************************************")
		
	# Perform analysis for Crab model
	pipeline("$CTOOLS/share/models/crab.xml")

	# Perform analysis for disk model
	pipeline("$CTOOLS/share/models/disk.xml")

	# Perform analysis for Gaussian model
	pipeline("$CTOOLS/share/models/gauss.xml")
	
	# Perform analysis for shell model
	pipeline("$CTOOLS/share/models/shell.xml")
