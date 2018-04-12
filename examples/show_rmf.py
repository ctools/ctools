#! /usr/bin/env python
# ==========================================================================
# Display Redistribution Matrix File
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
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()


# ======== #
# Plot RMF #
# ======== #
def plot_rmf(rmf, plotfile):
    """
    Plot Redistribution Matrix File

    Parameters
    ----------
    rmf : `~gammalib.GRmf`
        Redistribution Matrix File
    plotfile : str
        Plot filename
    """
    # Build title

    # Set energy range
    etrue_min = rmf.etrue().emin().log10TeV()
    etrue_max = rmf.etrue().emax().log10TeV()
    ereco_min = rmf.emeasured().emin().log10TeV()
    ereco_max = rmf.emeasured().emax().log10TeV()
    aspect    = (etrue_max-etrue_min)/(ereco_max-ereco_min)

    # Set number of bins
    netrue = rmf.etrue().size()
    nereco = rmf.emeasured().size()

    # Set image
    image = []
    for ireco in range(nereco-1,-1,-1):
        row = []
        for itrue in range(netrue):
            row.append(rmf[itrue,ireco])
        image.append(row)

    # Create figure
    fig = plt.figure()

    # Plot image
    c    = plt.imshow(image, extent=[etrue_min,etrue_max,ereco_min,ereco_max],
                      interpolation='nearest', aspect=aspect)
    cbar = plt.colorbar(c, orientation='vertical', shrink=0.8)

    # Plot title and axis
    plt.title('Redistribution Matrix File')
    plt.xlabel('log10(Etrue/TeV)')
    plt.ylabel('log10(Ereco/TeV)')

    # Show plots or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ======== #
# Show RMF #
# ======== #
def show_rmf():
    """
    Show Redistribution Matrix File
    """
    # Set usage string
    usage = 'show_rmf.py [-p plotfile] file'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Get RMF
    rmf = gammalib.GRmf(args[0])
    
    # Plot RMF
    plot_rmf(rmf, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show RMF
    show_rmf()
