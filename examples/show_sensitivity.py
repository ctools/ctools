#! /usr/bin/env python
# ==========================================================================
# This script displays a differential sensitivity plot generated by cssens.
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
import sys
import math
import gammalib
try:
    import matplotlib.pyplot as plt
except:
    sys.exit("This script needs matplotlib")


# ===================== #
# Read sensitivity data #
# ===================== #
def read_sensitivity(filename):
    """
    Read sensitivity information from CSV file.
    """
    # Read filename
    csv = gammalib.GCsv(filename,',')

    # Create dictionnary
    sensitivity = {}
    for column in range(csv.ncols()):
        name   = csv[0,column]
        values = []
        for row in range(csv.nrows()-1):
            values.append(float(csv[row+1,column]))
        sensitivity[name] = values

    # Check where file contains differential or integral sensitivity
    mode = "Integral"
    if sensitivity.has_key("emax"):
        emax_ref = -1.0
        for value in sensitivity["emax"]:
            if emax_ref < 0.0:
                emax_ref = value
            elif emax_ref != value:
                mode = "Differential"
                break

    # Add mode
    sensitivity["mode"] = mode

    # Add linear energy values
    if mode == "Differential":
        if sensitivity.has_key("loge"):
            name   = 'energy'
            values = []
            for value in sensitivity["loge"]:
                values.append(math.pow(10.0, value))
            sensitivity[name] = values
    else:
        if sensitivity.has_key("emin"):
            name   = 'energy'
            values = []
            for value in sensitivity["emin"]:
                values.append(value)
            sensitivity[name] = values

    # Return
    return sensitivity


# ===================== #
# Show sensitivity data #
# ===================== #
def show_sensitivity(sensitivity, filename):
    """
    Plot sensitivity.
    """
    # Build title
    title = sensitivity['mode'] + ' sensitivity ('+filename+')'

    # Create figure
    plt.figure()

    # Set plot attributes
    plt.title(title)
    plt.loglog()
    plt.grid()

    # Show differential sensitivity
    plt.plot(sensitivity['energy'], sensitivity['eflux'], 'ro-')

    # Set labels
    plt.xlabel("Energy (TeV)")
    if sensitivity['mode'] == "Differential":
        plt.ylabel(r"E$\times$ F(E) (erg cm$^{-2}$ s$^{-1}$)")
    else:
        plt.ylabel(r"E$\times$ F($>$E) (erg cm$^{-2}$ s$^{-1}$)")

    # Show plot
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    """
    # Print usage information
    usage = "Usage: show_sensitivity [filename]"
    if len(sys.argv) < 1:
        sys.stdout.write(usage+"\n")
        sys.exit()

    # Extract parameters
    if len(sys.argv) == 1:
        filename = 'sensitivity.dat'
    else:
        filename = sys.argv[1]

    # Read sensitivity data
    sensitivity = read_sensitivity(filename)

    # Show sensitivity data
    show_sensitivity(sensitivity, filename)
