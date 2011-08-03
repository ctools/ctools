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
def create_ts(loge, emin, emax, ntrials=3, duration=180000.0, \
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
	tsdist["irf"]      = "kb_E_50h_v3"
	tsdist["type"]     = "point"
	tsdist["index"]    = -2.48
	tsdist["offset"]   = 0.0
	tsdist["bkg"]      = "$GAMMALIB/share/models/bkg_kb_E_50h_v3.txt"
	tsdist["emin"]     = emin
	tsdist["emax"]     = emax
	tsdist["enumbins"] = enumbins
	tsdist["duration"] = duration
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
	usage = "ts_distributions [max_threads]"
	if len(sys.argv) < 1 or len(sys.argv) > 2:
		print usage
		sys.exit()

	# Set maximum number of threads (default: 1)
	if len(sys.argv) == 2:
		max_threads = int(sys.argv[1])
	else:
		max_threads = 1

	# Loop over energy bands. The energy bands are those that are also
	# used for sensitivity computation.
	for ieng in range(21):
		
		# Set energies
		loge  = -1.7 + ieng * 0.2
		emean = pow(10.0, loge)
		emin  = pow(10.0, loge-0.1)
		emax  = pow(10.0, loge+0.1)

		# Processing support?
		if has_processing:

			# Wait until one thread has finished
			while len(processing.activeChildren()) >= max_threads:
				time.sleep(60)

			# Set arguments
			args   = (loge, emin, emax)
			kwargs = {}

			# Generate pull distribution
			p = processing.Process(target=create_ts, args=args, kwargs=kwargs)
			p.start()
			print "Process emin=%.4f emax=%.4f started." % (emin, emax)

			# Wait a short time to allow process to start
			time.sleep(1)
		
		# ... no
		else:
			create_ts(loge, emin, emax)
	
	# Processing support
	if has_processing:
	
		# Wait until all threads finished
		while len(processing.activeChildren()) > 0:
			time.sleep(60)
