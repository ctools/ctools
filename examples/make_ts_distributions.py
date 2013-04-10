#! /usr/bin/env python
# ==========================================================================
# This script generates TS distributions for the standard energy bands that
# are also used for sensitivity computation. The script makes use of the
# "processing" Python module, if available, that can be used for parallel
# computing on multiple cores/CPUs.
#
# This script provides an illustration of how a ctools or cscript can be
# used from within Python. It further illustrates the possibility of
# parallel computing using the "processing" Python module.
#
# Required 3rd party modules:
# - processing (optional)
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
from cstsdist import *
import sys
import time

# Try importing processing module
has_processing = False
try:
	import processing
	has_processing = True
except:
	pass


# ====================== #
# Create TS distribution #
# ====================== #
def create_ts(loge, emin, emax, ntrials=100, duration=180000.0, \
			  enumbins=0, log=False):
	"""
	Create TS distribution.
	
	Parameters:
	 loge - Logarithm of mean energy in TeV
	 emin - Minimum energy (TeV)
	 emax - Maximum energy (TeV)
	Keywords:
	 log  - Create log file(s)
	"""
	# Generate output filename
	outfile = "ts_"+str(loge)+".dat"
	
	# Setup cstsdist tool
	tsdist = cstsdist()
	tsdist["outfile"]  = outfile
	tsdist["ntrials"]  = ntrials
	tsdist["caldb"]    = "$GAMMALIB/share/caldb/cta"
	tsdist["irf"]      = "cta_dummy_irf"
	tsdist["type"]     = "point"
	tsdist["index"]    = -2.48
	tsdist["offset"]   = 0.0
	tsdist["bkg"]      = "$GAMMALIB/share/models/bkg_dummy.txt"
	tsdist["emin"]     = float(emin)
	tsdist["emax"]     = float(emax)
	tsdist["enumbins"] = int(enumbins)
	tsdist["duration"] = float(duration)
	tsdist["rad"]      = 5.0
	tsdist["npix"]     = 200
	tsdist["binsz"]    = 0.05
		
	# Optionally open the log file
	if log:
		tsdist.logFileOpen()
	
	# Run tool
	tsdist.run()
	
	# Return
	return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
	"""
	Create TS distribution in a number of energy bands.
	"""
	# Get input arguments
	usage = "make_ts_distributions [-n ntrials] [-e enumbins] [-m max_threads]"
	if len(sys.argv) < 1:
		print(usage)
		sys.exit()
	
	# Set default parameters
	ntrials     = 100
	enumbins    = 0
	duration    = 180000.0
	max_threads = 1
	
	# Parameter dictionnary
	pars = [{'option': '-n', 'value': ntrials}, \
	        {'option': '-e', 'value': enumbins}, \
	        {'option': '-d', 'value': duration}, \
			{'option': '-m', 'value': max_threads}]
	
	# Gather parameters from command line
	i = 1
	while i < len(sys.argv):
		
		# Search for options
		for par in pars:
			if sys.argv[i] == par['option']:
				if len(sys.argv) > i+1:
					i      += 1
					try:
						par['value'] = int(sys.argv[i])
					except:
						print(usage)
						sys.exit()
				else:
					print(usage)
					sys.exit()
		
		# Next item
		i += 1
	
	# Recover parameters
	ntrials     = pars[0]['value']
	enumbins    = pars[1]['value']
	duration    = pars[2]['value']
	max_threads = pars[3]['value']
	#print ntrials
	#print enumbins
	#print duration
	#print max_threads
	
	# Loop over energy bands. The energy bands are those that are also
	# used for sensitivity computation.
	for ieng in range(21):
		
		# Set energies
		loge  = -1.7 + ieng * 0.2
		emean = pow(10.0, loge)
		emin  = pow(10.0, loge-0.1)
		emax  = pow(10.0, loge+0.1)
		if loge < 0:
			loge = "m"+str(abs(loge))
		else:
			loge = "p"+str(abs(loge))

		# Processing support?
		if has_processing:

			# Wait until one thread has finished
			while len(processing.activeChildren()) >= max_threads:
				time.sleep(10)

			# Set arguments
			args   = (loge, emin, emax)
			kwargs = {'ntrials': ntrials, 'enumbins': enumbins, \
			          'duration': duration}

			# Generate pull distribution
			p = processing.Process(target=create_ts, args=args, kwargs=kwargs)
			p.start()
			print("Process emin=%.4f emax=%.4f started." % (emin, emax))

			# Wait a short time to allow process to start
			time.sleep(1)
		
		# ... no
		else:
			create_ts(loge, emin, emax, ntrials=ntrials, enumbins=enumbins, \
                      duration=duration)
	
	# Processing support
	if has_processing:
	
		# Wait until all threads finished
		while len(processing.activeChildren()) > 0:
			time.sleep(10)
