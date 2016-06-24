#! /usr/bin/env python
# ==========================================================================
# Shows the pull histogram
#
# Copyright (C) 2011-2016 Juergen Knoedlseder
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
try:
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
except:
    print('This script needs the "matplotlib" module')
    sys.exit()
try:
    import numpy as np
except:
    print('This script needs the "numpy" module')
    sys.exit()


# ====================== #
# Read pull distribution #
# ====================== #
def read_pull(filename, parname):
    """
    Read pull distribution
    
    Parameters
    ----------
    filename : str
        Pull distribution ASCII file
    parname : str
        Parameter
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
                print('ERROR: Parameter "'+parname+'" not found in list:')
                for p in row:
                    print('       "'+p+'"')
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
    usage = 'Usage: show_pull_histogram filename parname [bins] [file]'
    if len(sys.argv) < 3:
        print(usage)
        sys.exit()

    # Extract parameters
    filename = sys.argv[1]
    parname  = sys.argv[2]
    if len(sys.argv) > 3:
        nbins = int(sys.argv[3])
    else:
        nbins = 50
    if len(sys.argv) == 5:
        plotfile = sys.argv[4]
    else:
        plotfile = ''

    # Read values from CSV file
    values = read_pull(filename, parname)

    # Create histogram
    n, bins, patches = plt.hist(values, nbins, range=[-4.0,4.0],
                                normed=True, facecolor='green')

    # Create expected distribution
    y = mlab.normpdf(bins, 0.0, 1.0)
    l = plt.plot(bins, y, 'r-', linewidth=2)

    # Set parname for plotting
    name = parname

    # Set plot
    plt.xlabel('Pull ('+name+')')
    plt.ylabel('Arbitrary units')
    plt.title(name)
    plt.grid(True)

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()
