#! /usr/bin/env python
# ==========================================================================
# This script shows how to plot a butterfly created with ctbutterfly 
#
# Copyright (C) 2014-2015 Michael Mayer
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

import gammalib
import sys
try:
    import matplotlib.pyplot as plt
except:
    sys.exit("This script needs matplotlib")

# =============== #
# Script entry    #
# =============== #    
if __name__ == "__main__":

    # check for given butterfly file
    if not len(sys.argv) == 2:
        sys.exit("Usage: show_butterfly.py butterfly.txt")

    # Read given butterfly file    
    filename = sys.argv[1]
    csv = gammalib.GCsv(filename)

    # initialise arrays to be filled
    butterfly_x = []
    butterfly_y = []

    line_x = []
    line_y = []

    # loop over rows of the file
    nrows = csv.nrows()
    for row in range(nrows):
        # Compute upper edge of confidence band
        butterfly_x.append(csv.real(row,0))
        butterfly_y.append(csv.real(row,1)+csv.real(row,2))

        # Set line values
        line_x.append(csv.real(row,0))
        line_y.append(csv.real(row,1))

    # Loop over the rows backwards to compute the lower edge
    # of the confidence band    
    for row in range(nrows):
        index = nrows - 1 - row

        butterfly_x.append(csv.real(index,0))
        low_error = csv.real(index,1)-csv.real(index,2)
        if low_error < 1e-26:
            low_error = 1e-26
        butterfly_y.append(low_error)   
    
    # plot the butterfly and spectral line       
    plt.figure()
    plt.ylim([1e-26,1e-14])
    plt.loglog()
    plt.grid()
    plt.plot(line_x,line_y,color='black',ls='-')
    plt.fill(butterfly_x,butterfly_y,color='green',alpha=0.5)

    plt.xlabel("Energy [MeV]")
    plt.ylabel(r"dN/dE [MeV$^{-1}$ s$^{-1}$ cm$^{-2}$]")    
    plt.show()
