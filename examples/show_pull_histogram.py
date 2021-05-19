#! /usr/bin/env python
# ==========================================================================
# Shows the pull histogram
#
# Copyright (C) 2011-2021 Juergen Knoedlseder
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


# ====================== #
# Read pull distribution #
# ====================== #
def read_pull_distribution(filename, parname):
    """
    Read pull distribution

    Parameters
    ----------
    filename : str
        Name of FITS or CSV file
    parname : str
        Parameter name

    Returns
    -------
    pulls : list of floats
        List with pull distribution
    """
    # Read pull distribution dependent on file type
    fname = gammalib.GFilename(filename)
    if fname.is_fits():
        pulls = read_pull_distribution_fits(filename, parname)
    else:
        pulls = cscripts.ioutils.read_pull_values(filename, parname)

    # Return pulls
    return pulls


# ===================================== #
# Read pull distribution from FITS file #
# ===================================== #
def read_pull_distribution_fits(filename, parname):
    """
    Read pull distribution from FITS file

    Parameters
    ----------
    filename : str
        Name of FITS file
    parname : str
        Parameter name

    Returns
    -------
    pulls : list of floats
        List with pull distribution
    """
    # Open FITS file
    fits = gammalib.GFits(filename)

    # Get pull distribution table
    table = fits.table('PULL_DISTRIBUTION')

    # Get relevant column
    column = table[parname]

    # Initialise vector
    pulls = []

    # Fill vectors
    nrows = table.nrows()
    for row in range(nrows):
        pulls.append(column[row])

    # Return
    return pulls


# =================== #
# Plot pull histogram #
# =================== #
def plot_pull_histogram(filename, parname, nbins, plotfile):
    """
    Plot pull histogram

    Parameters
    ----------
    filename : str
        Pull filename
    parname : str
        Parameter name
    nbins : int
        Number of hisogram bins
    plotfile : str
        Plot filename
    """
    # Read pull distribution from file
    pulls = read_pull_distribution(filename, parname)
    
    # Create Numpy array
    values = np.array(pulls)

    # Create histogram
    _, bins, _ = plt.hist(values, nbins, range=[-4.0,4.0],
                          normed=True, facecolor='green')

    # Create expected distribution
    y = np.exp(-0.5*bins*bins) / np.sqrt(2.0*np.pi)
    plt.plot(bins, y, 'r-', linewidth=2)

    # Set plot
    plt.xlabel('Pull ('+parname+')')
    plt.ylabel('Arbitrary units')
    plt.title(parname)
    plt.grid(True)

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# =================== #
# Show pull histogram #
# =================== #
def show_pull_histogram():
    """
    Show pull histogram
    """
    # Set usage string
    usage = 'show_pull_histogram.py [-n bins] [-p plotfile] file parameter'

    # Set default options
    options = [{'option': '-n', 'value': '50'},
               {'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    nbins    = int(options[0]['value'])
    plotfile = options[1]['value']

    # Plot pull histogram
    plot_pull_histogram(args[0], args[1], nbins, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show pull histogram
    show_pull_histogram()
