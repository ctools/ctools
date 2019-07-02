#! /usr/bin/env python
# ==========================================================================
# Shows butterfly diagram created with ctbutterfly
#
# Copyright (C) 2014-2017 Michael Mayer
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


def get_butterfly_file(filename):
    """
    Get butterfly diagram values

    Parameters
    ----------
    filename : str
        Butterfly CSV file
    
    Returns
    -------
    See `get_butterfly_csv`
    """
    # Read butterfly file
    csv = gammalib.GCsv(filename)

    return get_butterfly_csv(csv)


def get_butterfly_csv(csv):
    """
    Get butterfly diagram values from GCsv object

    Parameters
    ----------
    csv : `~gammalib.GCsv`
        GCsv object containing butterfly plot information

    Returns
    -------
    'dict' defining butterfly plot and best fit spectrum
    """
    # Initialise arrays to be filled
    btrfly = {
        'butterfly_x' : [],
        'butterfly_y' : [],
        'line_x'      : [],
        'line_y'      : []
    }

    # Loop over rows of the file
    nrows = csv.nrows()
    for row in range(nrows):

        # Get conversion coefficient
        conv = csv.real(row,0) * csv.real(row,0) * gammalib.MeV2erg

        # Compute upper edge of confidence band
        btrfly['butterfly_x'].append(csv.real(row,0)/1.0e6) # TeV
        btrfly['butterfly_y'].append(csv.real(row,2)*conv)

        # Set line values
        btrfly['line_x'].append(csv.real(row,0)/1.0e6) # TeV
        btrfly['line_y'].append(csv.real(row,1)*conv)

    # Loop over the rows backwards to compute the lower edge of the
    # confidence band
    for row in range(nrows):
        index = nrows - 1 - row
        conv  = csv.real(index,0) * csv.real(index,0) * gammalib.MeV2erg
        btrfly['butterfly_x'].append(csv.real(index,0)/1.0e6)
        low_error = max(csv.real(index,3)*conv, 1e-26)
        btrfly['butterfly_y'].append(low_error)

    return btrfly


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
    btrfly = get_butterfly_file(filename)
    
    # Plot the butterfly and spectral line       
    plt.figure()
    plt.loglog()
    plt.grid()
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)')
    
    # Plot the butterfly and spectral line
    plt.plot(btrfly['line_x'],btrfly['line_y'],color='black',ls='-')
    plt.fill(btrfly['butterfly_x'],btrfly['butterfly_y'],color='green',alpha=0.5)

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
    plot_butterfly(filename=args[0], plotfile=plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show butterfly diagram
    show_butterfly()
