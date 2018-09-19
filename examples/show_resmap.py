#! /usr/bin/env python
# ==========================================================================
# Shows histogram of residual map
#
# Copyright (C) 2018 Juergen Knoedlseder
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
import gammalib
import cscripts
try:
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
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


# ================= #
# Plot residual map #
# ================= #
def plot_resmap(filename, nbins, plotfile):
    """
    Plot pull histogram

    Parameters
    ----------
    filename : str
        Residual map filename
    nbins : int
        Number of hisogram bins
    plotfile : str
        Plot filename
    """
    # Open residual map
    map = gammalib.GSkyMap(filename)

    # Create array of values
    values = map.array()

    # Create histogram
    _, bins, _ = plt.hist(values, nbins, range=[-4.0,4.0],
                          normed=True, facecolor='green')

    # Create expected distribution
    y = mlab.normpdf(bins, 0.0, 1.0)
    plt.plot(bins, y, 'r-', linewidth=2)

    # Set plot
    plt.xlabel('Pull')
    plt.ylabel('Arbitrary units')
    #plt.title(parname)
    plt.grid(True)

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ================= #
# Show residual map #
# ================= #
def show_resmap():
    """
    Show residual map
    """
    # Set usage string
    usage = 'show_resmap.py [-n bins] [-p plotfile] file'

    # Set default options
    options = [{'option': '-n', 'value': '50'},
               {'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    nbins    = int(options[0]['value'])
    plotfile = options[1]['value']

    # Plot residual map
    plot_resmap(args[0], nbins, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show residual map
    show_resmap()
