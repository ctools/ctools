#! /usr/bin/env python
# ==========================================================================
# Display TS distribution generated by cstsdist
#
# Copyright (C) 2011-2021 Jurgen Knodlseder
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
import sys
import csv
import math
try:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()
try:
    import numpy as np
except ImportError:
    print('This script needs the "numpy" module')
    sys.exit()
import gammalib
import cscripts


# ==================== #
# Read TS distribution #
# ==================== #
def read_ts_distribution(filename):
    """
    Read TS distribution

    Parameters
    ----------
    filename : str
        Name of FITS or CSV file
    parname : str
        Parameter name

    Returns
    -------
    ts : list of floats
        List with TS distribution
    """
    # Read TS distribution dependent on file type
    fname = gammalib.GFilename(filename)
    if fname.is_fits():
        ts = read_ts_distribution_fits(filename)
    else:
        ts = read_ts_distribution_csv(filename)

    # Return ts
    return ts


# =================================== #
# Read TS distribution from FITS file #
# =================================== #
def read_ts_distribution_fits(filename):
    """
    Read TS distribution from FITS file

    Parameters
    ----------
    filename : str
        Name of FITS file
    parname : str
        Parameter name

    Returns
    -------
    ts : list of floats
        List with TS distribution
    """
    # Open FITS file
    fits = gammalib.GFits(filename)

    # Get TS distribution table
    table = fits.table('TS_DISTRIBUTION')

    # Get relevant column
    column = table['TS']

    # Initialise vector
    ts = []

    # Fill vectors
    nrows = table.nrows()
    for row in range(nrows):
        ts.append(column[row])

    # Return TS distribution
    return ts


# ==================== #
# Read TS distribution #
# ==================== #
def read_ts_distribution_csv(filename, tsname='TS'):
    """
    Read TS distribution values from CSV file

    Parameters
    ----------
    filename : str
        Name of CSV file
    tsname : str, optional
        Column name of TS

    Returns
    -------
    values : list of floats
        List with TS distribution
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
                print('ERROR: Column "'+tsname+'" not found in file')
                print(row)
                raise NameError(tsname)

        # Handle data rows
        else:
            values.append(float(row[index]))

        # Flag that first row has been passed
        first = False

    # Return TS values
    return values


# =========== #
# Compute erf #
# =========== #
def erf(x):
    """
    Compute error function

    Parameters
    ----------
    x : float
        Value

    Returns
    -------
    erf(-x) : float
        Error function
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
def dist_cdf(values, nbins, title, plotfile, expect=''):
    """
    Display cumulative TS distribution and overplots expectation

    Parameters
    ----------
    values : float
        TS values
    nbins : int
        Number of bins
    title : str
        Plot title
    plotfile : str
        Plot file name
    expect : str, optional
        If 'show', show expectation
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
        #x_val = 0.5*(x_min+x_max)
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
    plt.semilogy(hist_x, hist_y,      'r-', linewidth=2, label='simulations')
    if expect == 'show':
        plt.semilogy(hist_x, hist_expect, 'k-', linewidth=2, label='expectation')
        plt.semilogy(hist_x, hist_y,      'r-', linewidth=2)

    # Set plot
    plt.xlabel('TS')
    plt.ylabel('Number of trials > TS')
    plt.title(title)
    plt.grid(True)
    plt.legend(loc='upper right')

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ======================================== #
# Display probability density distribution #
# ======================================== #
def dist_pdf(values, nbins, title, plotfile, expect=''):
    """
    Display probability density TS distribution and overplots expectation

    Parameters
    ----------
    values : float
        TS values
    nbins : int
        Number of bins
    title : str
        Plot title
    plotfile : str
        Plot file name
    expect : str, optional
        If 'show', show expectation
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
    _, bins, _ = plt.hist(values, nbins, range=[min_value, max_value],
                          align='mid', facecolor='green',
                          log=True, label='simulations')

    # Create expected distribution (only for positive TS). We compute the
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
    for i in bins:
        ts = i + 0.5*width
        if ts > 0:
            y  = norm*math.pow(ts,-0.5)*math.exp(-0.5*ts)
            x.append(ts)
            y_full.append(y)
            y_half.append(0.5*y)

    # Show expected distribution (half)
    if expect == 'show':
        plt.semilogy(x, y_half, 'ro', linewidth=2, label='expectation')

    # Set plot
    plt.xlabel('TS')
    plt.ylabel('Number of trials')
    plt.title(title)
    plt.grid(True)
    plt.legend(loc='upper right')

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ==================== #
# Show TS distribution #
# ==================== #
def show_ts_distribution():
    """
    Show TS distribution
    """
    # Set usage string
    usage = ('show_ts_distribution.py [-n bins] [-c column] [-title title] '
             '[-e expect] [-t type] [-p plotfile] file')

    # Set default options
    options = [{'option': '-n', 'value': 30},
               {'option': '-c', 'value': 'TS'},
               {'option': '-title', 'value': 'TS distribution'},
               {'option': '-e', 'value': ''},
               {'option': '-t', 'value': 'pdf'},
               {'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    nbins    = int(options[0]['value'])
    tsname   = options[1]['value']
    title    = options[2]['value']
    expect   = options[3]['value']
    disttype = options[4]['value']
    plotfile = options[5]['value']

    # Read TS distribution from file
    ts = read_ts_distribution(args[0])

    # Create Numpy array
    values = np.array(ts)

    # Show histogram
    if disttype == 'pdf':
        dist_pdf(values, nbins, title, plotfile, expect=expect)
    else:
        dist_cdf(values, nbins, title, plotfile, expect=expect)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show TS distribution
    show_ts_distribution()
