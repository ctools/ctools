#! /usr/bin/env python
# ==========================================================================
# Shows the evolution of the mean and rms pull
#
# Copyright (C) 2015-2016 Juergen Knoedlseder
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
import cscripts
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


# =================== #
# Plot pull evolution #
# =================== #
def plot_pull_evolution(filename, parname, plotfile):
    """
    Plot pull evolution

    Parameters
    ----------
    filename : str
        Pull filename
    parname : str
        Parameter name
    plotfile : str
        Plot filename
    """
    # Read values from CSV file
    values = np.array(cscripts.ioutils.read_pull_values(filename, parname))

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
    plt.semilogx(samples, mean, 'r-', label='Mean')
    plt.semilogx(samples, rms, 'b-', label='Std. deviation')

    # Set plot
    plt.xlabel('Samples')
    plt.ylabel('Value')
    plt.title(parname)
    plt.grid(True)

    # Set legend
    plt.legend(loc='upper left')

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# =================== #
# Show pull evolution #
# =================== #
def show_pull_evolution():
    """
    Show pull evolution
    """
    # Set usage string
    usage = 'show_pull_evolution.py [-p plotfile] file parameter'

    # Set default options
    options = [{'option': '-p',   'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Plot pull evolution
    plot_pull_evolution(args[0], args[1], plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show pull evolution
    show_pull_evolution()
