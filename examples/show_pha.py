#! /usr/bin/env python
# ==========================================================================
# Show a spectrum in PHA format
#
# Copyright (C) 2013-2016 Juergen Knoedlseder
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
import math
import gammalib
import cscripts
try:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()


# ================= #
# Plot PHA spectrum #
# ================= #
def plot_pha(pha, plotfile):
    """
    Plot PHA spectrum

    Parameters
    ----------
    spectrum : `~gammalib.GPha`
        PHA file
    plotfile : str
        Plot file name
    """
    # Create figure
    plt.figure(1)
    plt.title('PHA spectrum (' + pha.filename() + ')')

    # Generate energy vector
    energy   = []
    channels = False
    if pha.ebounds().size() > 0:
        for i in range(pha.ebounds().size()):
            energy.append(pha.ebounds().elogmean(i).TeV())
    else:
        for i in range(pha.size()):
            energy.append(float(i+1))
            channels = True

    # Generate data and error vector
    counts = [0.0 for i in range(pha.size())]
    for i in range(pha.size()):
        counts[i] = pha[i]
    error = [math.sqrt(c) for c in counts]

    # Plot data
    if channels:
        plt.semilogy(energy, counts, 'ro')
    else:
        plt.loglog(energy, counts, 'ro')
    plt.errorbar(energy, counts, error, fmt='ro', ecolor='r')

    # Set axes
    if channels:
        plt.xlabel('Channels')
    else:
        plt.xlabel('Energy (TeV)')
    plt.ylabel('Counts')

    # Show spectrum or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ================= #
# Show PHA spectrum #
# ================= #
def show_pha():
    """
    Show PHA spectrum
    """
    # Set usage string
    usage = 'show_pha.py [-p plotfile] file'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Load PHA spectrum
    pha = gammalib.GPha(args[0])

    # Plot PHA spectrum
    plot_pha(pha, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show PHA spectrum
    show_pha()
