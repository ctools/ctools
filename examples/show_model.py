#! /usr/bin/env python
# ==========================================================================
# This script shows a model 
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
import gammalib
import sys
import math
try:
    import matplotlib.pyplot as plt
except:
    sys.exit("This script needs matplotlib")


# ======================= #
# Plot spectral component #
# ======================= #    
def plot_spectrum(model, emin=0.01, emax=100.0, enumbins=100):
    """
    """
    # Get spectral component
    spectrum = model.spectral()

    # Setup energy axis
    e_min   = gammalib.GEnergy(emin, "TeV")
    e_max   = gammalib.GEnergy(emax, "TeV")
    ebounds = gammalib.GEbounds(enumbins, e_min, e_max)

    # Setup model plot
    x   = []
    y   = []
    min = 1.0e30
    max = 0.0
    for i in range(enumbins):
        energy = ebounds.elogmean(i)
        value  = spectrum.eval(energy, gammalib.GTime())
        if value > max:
            max = value
        if value < min:
            min = value
        x.append(energy.TeV())
        y.append(value)

    # Set minimum and maximum
    max = math.pow(10.0, math.log10(max)+1)
    if min < 0.0001*max:
        min = 0.0001*max
    min = math.pow(10.0, math.log10(min)-1)


    # Show spectrum
    plt.figure()
    plt.title(model.name()+" ("+spectrum.type()+")")
    plt.loglog()
    plt.grid()
    plt.loglog(x, y, color='red')
    plt.xlabel("Energy [TeV]")
    plt.ylabel(r"dN/dE [ph s$^{-1}$ cm$^{-2}$ MeV$^{-1}$]")
    plt.ylim([min,max]) 
    plt.show()

    # Return
    return


# ============= #
# Script entry  #
# ============= #    
if __name__ == "__main__":

    # Check for model filename
    if len(sys.argv) < 2:
        sys.exit("Usage: show_model.py model.xml [index]")

    # Get optional model index
    if len(sys.argv) == 3:
        index = int(sys.argv[2])
    else:
        index = 0

    # Read models XML file
    filename = sys.argv[1]
    models   = gammalib.GModels(filename)

    # Extract relevant model
    model = models[index]

    # Plot spectrum
    plot_spectrum(model)
