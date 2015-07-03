#! /usr/bin/env python
# ==========================================================================
# This script displays the TS distribution values generated using cstsdist.
#
# Required 3rd party modules:
# - matplotlib
# - numpy
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
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv
import math


# ==================== #
# Read TS distribution #
# ==================== #
def read_ts(filename, tsname="TS"):
	"""
	Read TS distribution values from CSV file.
	"""
	# Initialise list
	values = []
	
	# Open reader
	reader = csv.reader(open(filename, 'r'), delimiter=',')
	
	# Read rows
	first = True
	index = -1
	for row in reader:
		
		# Get column index if first row
		if first:
			try:
				index = row.index(tsname)
			except:
				sys.stdout.write('ERROR: Column "TS" not found in file\n')
				sys.stdout.write(row+"\n")
				raise NameError("TS")

		# Handle data rows
		else:
			values.append(float(row[index]))
		
		# Flag that first row has been passed
		first = False
	
	# Create numpy array
	a = np.array(values)
	
	# Return array
	return a


# =========== #
# Compute erf #
# =========== #
def erf(x):
	"""
	Compute error function.
	"""
	# save the sign of x
	sign = 1
	if x < 0: 
		sign = -1
	x = abs(x)

	# constants
	a1 =  0.254829592
	a2 = -0.284496736
	a3 =  1.421413741
	a4 = -1.453152027
	a5 =  1.061405429
	p  =  0.3275911
	
	# A&S formula 7.1.26
	t = 1.0/(1.0 + p*x)
	y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
	
	# Return result
	return sign*y # erf(-x) = -erf(x)


# =============================== #
# Display cumulative distribution #
# =============================== #
def dist_cdf(values, nbins, title):
	"""
	Displays cumulative TS distribution and overplots expectation.
	"""
	# Set range and adapt number of bins. We make sure that we have a
	# histogram bin centred on 0 that will capture the TS=0 part.
	max_value = max(values)
	binsize   = max_value/nbins
	min_value = 0.5*binsize
	while min_value > min(values):
		min_value -= binsize
		nbins     += 1
	
	# Create histogram
	hist_x      = []
	hist_y      = []
	hist_expect = []
	norm        = 0.5*len(values)
	for i in range(nbins):
		x_min = min_value+i*binsize
		x_max = min_value+(i+1)*binsize
		x_val = 0.5*(x_min+x_max)
		y_val = 0.0
		for ts in values:
			if ts > x_min:
				y_val += 1.0
		hist_x.append(x_min)
		hist_y.append(y_val)
		hist_x.append(x_max)
		hist_y.append(y_val)
		if x_min > 0.0:
			value = norm*erf(math.sqrt(0.5*x_min)) + norm
			hist_expect.append(len(values)-value)
			hist_expect.append(len(values)-value)
		else:
			hist_expect.append(norm)
			hist_expect.append(norm)

	# Show expected distribution (half)
	l = plt.semilogy(hist_x, hist_y,      'r-', linewidth=2, label="simulations")
	l = plt.semilogy(hist_x, hist_expect, 'k-', linewidth=2, label="expectation")
	l = plt.semilogy(hist_x, hist_y,      'r-', linewidth=2)

	# Set plot
	plt.xlabel('TS')
	plt.ylabel('Number of trials > TS')
	plt.title(title)
	plt.grid(True)
	plt.legend(loc="upper right")

	# Show histogram
	plt.show()
	
	# Return
	return


# ======================================== #
# Display probability density distribution #
# ======================================== #
def dist_pdf(values, nbins, title):
	"""
	Display probability density TS distribution and overplots
	expectation.
	"""
	# Set range and adapt number of bins. We make sure that we have a
	# histogram bin centred on 0 that will capture the TS=0 part.
	max_value = max(values)
	binsize   = max_value/nbins
	min_value = 0.5*binsize
	while min_value > min(values):
		min_value -= binsize
		nbins     += 1
	
	# Create histogram
	n, bins, patches = plt.hist(values, nbins, range=[min_value, max_value],
	                            align='mid', facecolor='green',
								log=True, label="simulations")

	# Create expected distribution (onyl for positive TS). We compute the
	# full and half of the distribution, as for positively constrained
	# quantities half of the values have TS=0
	x      = []
	y_full = []
	y_half = []
	width  = bins[1]-bins[0]
	norm   = len(values)/(math.sqrt(2.0)*math.sqrt(math.pi))*width
	x.append(0.0)
	y_full.append(len(values))
	y_half.append(0.5*len(values))
	for bin in bins:
		ts = bin + 0.5*width
		if ts > 0:
			y  = norm*math.pow(ts,-0.5)*math.exp(-0.5*ts)
			x.append(ts)
			y_full.append(y)
			y_half.append(0.5*y)
	
	# Show expected distribution (half)
	l = plt.semilogy(x, y_half, 'ro', linewidth=2, label="expectation")

	# Set plot
	plt.xlabel('TS')
	plt.ylabel('Number of trials')
	plt.title(title)
	plt.grid(True)
	plt.legend(loc="upper right")

	# Show histogram
	plt.show()
		
	
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
	"""
	Display TS distribution generated using cstsdist using matplotlib.
	"""
	# Print usage information
	usage = "Usage: show_ts_distribution filename" \
	        " [-n bins] [-c column] [-t title] [-p plot]"
	if len(sys.argv) < 2:
		sys.stdout.write(usage+"\n")
		sys.exit()

	# Extract parameters
	filename = sys.argv[1]

	# Set default parameters
	nbins  = 30
	tsname = "TS"
	title  = "TS distribution"
	plot   = "pdf"

	# Parameter dictionnary
	pars = [{'option': '-n', 'value': nbins},
	        {'option': '-c', 'value': tsname},
	        {'option': '-t', 'value': title},
			{'option': '-p', 'value': plot}]
	
	# Gather parameters from command line
	i = 2
	while i < len(sys.argv):
		
		# Search for options
		for par in pars:
			if sys.argv[i] == par['option']:
				if len(sys.argv) > i+1:
					i      += 1
					try:
						par['value'] = sys.argv[i]
					except:
						sys.stdout.write(usage+"\n")
						sys.exit()
				else:
					sys.stdout.write(usage+"\n")
					sys.exit()
		
		# Next item
		i += 1
	
	# Recover parameters
	nbins  = int(pars[0]['value'])
	tsname = pars[1]['value']
	title  = pars[2]['value']
	plot   = pars[3]['value']
	
	# Read values from CSV file
	values = read_ts(filename, tsname=tsname)
	sys.stdout.write(len(values), "values read.\n")
	
	# Show histogram
	if plot == "pdf":
		dist_pdf(values, nbins, title)
	else:
		dist_cdf(values, nbins, title)
