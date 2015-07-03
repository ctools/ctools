#! /usr/bin/env python
# ==========================================================================
# This script displays a spectrum in PHA format.
#
# Required 3rd party modules:
# - matplotlib
# - numpy
#
# Copyright (C) 2013 Juergen Knoedlseder
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


# ================= #
# Show PHA spectrum #
# ================= #
def show_pha(pha):
    """
    """
    # Create figure
    plt.figure(1)
    plt.title("PHA spectrum (" + pha.filename() + ")")

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
    plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

    # Set axes
    if channels:
        plt.xlabel("Channels")
    else:
        plt.xlabel("Energy (TeV)")
    plt.ylabel("Counts")

    # Show plot
    plt.show()
    
    # Return
    #return
    
    
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Display spectrum in PHA format.
    """
    # Print usage information
    usage = "Usage: show_pha filename" \
            " [-n bins] [-c column] [-t title] [-p plot]"
    if len(sys.argv) < 2:
        print(usage)
        sys.exit()

    # Extract parameters
    filename = sys.argv[1]

    # Load PHA spectrum
    pha = gammalib.GPha(filename)
    #print(pha)

    # Try importing matplotlib
    try:
        import matplotlib.pyplot as plt
        has_matplotlib = True
    except ImportError:
        print("Matplotlib is not (correctly) installed on your system.")
        has_matplotlib = False

    # Show PHA spectrum
    if has_matplotlib:
        show_pha(pha)
    