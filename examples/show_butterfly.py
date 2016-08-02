#! /usr/bin/env python
# ==========================================================================
# Shows butterfly diagram created with ctbutterfly
#
# Copyright (C) 2014-2016 Michael Mayer
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
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()


# ============== #
# Plot butterfly #
# ============== #
def plot_butterfly(filename, plotfile):
    """
    Plot butterfly diagram

    Parameters
    ----------
    filename : str
        Butterfly CSV file
    plotfile : str
        Plot filename
    """
    # Read butterfly file
    csv = gammalib.GCsv(filename)

    # Initialise arrays to be filled
    butterfly_x = []
    butterfly_y = []
    line_x      = []
    line_y      = []

    # Loop over rows of the file
    nrows = csv.nrows()
    for row in range(nrows):

        # Compute upper edge of confidence band
        butterfly_x.append(csv.real(row,0))
        butterfly_y.append(csv.real(row,2))

        # Set line values
        line_x.append(csv.real(row,0))
        line_y.append(csv.real(row,1))

    # Loop over the rows backwards to compute the lower edge of the
    # confidence band
    for row in range(nrows):
        index = nrows - 1 - row
        butterfly_x.append(csv.real(index,0))
        low_error = max(csv.real(index,3), 1e-26)
        butterfly_y.append(low_error)
    
    # Plot the butterfly and spectral line       
    plt.figure()
    plt.ylim([1e-26,1e-14])
    plt.loglog()
    plt.grid()
    plt.plot(line_x,line_y,color='black',ls='-')
    plt.fill(butterfly_x,butterfly_y,color='green',alpha=0.5)
    plt.xlabel('Energy (MeV)')
    plt.ylabel(r'E $\times$ dN/dE (MeV$^{-1}$ s$^{-1}$ cm$^{-2}$)')

    # Show spectrum or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ====================== #
# Show butterfly diagram #
# ====================== #
def show_butterfly():
    """
    Show butterfly diagram
    """
    # Set usage string
    usage = 'show_butterfly.py [-p plotfile] file'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Plot butterfly diagram
    plot_butterfly(args[0], plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show butterfly diagram
    show_butterfly()
