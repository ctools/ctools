#! /usr/bin/env python
# ==========================================================================
# This script displays the evolution of the pull distributions as function
# of the number of samples. The pull distributions are generated using the
# cspull script. The pull is defined by
#
#           (fitted value - simulated value) / estimated error
#
# and the pull histogram allows to assess possible biases in the fit
# parameters and errors.
#
# Required 3rd party modules:
# - matplotlib
# - numpy
#
# Copyright (C) 2015 Juergen Knoedlseder
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
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv


# ====================== #
# Read pull distribution #
# ====================== #
def read_pull(filename, parname):
    """
    Read pull distribution.
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
                index = row.index(parname)
            except:
                sys.stdout.write('ERROR: Parameter "'+parname+'" not found in list:\n')
                for p in row:
                    sys.stdout.write('       "'+p+'"\n')
                raise NameError(parname)

        # Handle data rows
        else:
            values.append(float(row[index]))

        # Flag that first row has been passed
        first = False

    # Create numpy array
    a = np.array(values)

    # Return array
    return a


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Print usage information
    usage = "Usage: show_pull_evolution filename parname"
    if len(sys.argv) < 3:
        sys.stdout.write(usage+"\n")
        sys.exit()

    # Extract parameters
    filename = sys.argv[1]
    parname  = sys.argv[2]

    # Read values from CSV file
    values = read_pull(filename, parname)

    # Compute mean and rms
    nsamples = len(values)
    samples  = []
    mean     = []
    rms      = []
    for i in range(1,nsamples):
        samples.append(float(i+1))
        mean.append(np.mean(values[0:i]))
        rms.append(np.std(values[0:i]))

    # Plot mean and rms evolution
    plt.semilogx(samples, mean, 'r-', label="Mean")
    plt.semilogx(samples, rms, 'b-', label="Std. deviation")

    # Set plot
    plt.xlabel("Samples")
    plt.ylabel("Value")
    plt.title(parname)
    plt.grid(True)

    # Set legend
    plt.legend(loc="upper left")

    # Show histogram
    plt.show()
