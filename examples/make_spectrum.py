#! /usr/bin/env python
# ==========================================================================
# This script shows how to make a spectrum using obsutils.
#
# Copyright (C) 2014 Juergen Knoedlseder
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
import ctools
from ctools import obsutils
try:
	import matplotlib.pyplot as plt
	has_matplotlib = True
except:
	has_matplotlib = False


# =============== #
# Make a spectrum #
# =============== #
def make_spectrum():
    """
    """
    # Set script parameters
    model_name  = "${CTOOLS}/share/models/crab.xml"
    caldb       = "prod2"
    irf         = "South_50h"
    ra          =   83.63
    dec         =   22.01
    rad_sim     =    3.0
    tstart      =    0.0
    tstop       = 1800.0
    emin        =    0.1
    emax        =  100.0

    # Simulate events
    sim = ctools.ctobssim()
    sim["inmodel"].filename(model_name)
    sim["caldb"].string(caldb)
    sim["irf"].string(irf)
    sim["ra"].real(ra)
    sim["dec"].real(dec)
    sim["rad"].real(rad_sim)
    sim["tmin"].real(tstart)
    sim["tmax"].real(tstop)
    sim["emin"].real(emin)
    sim["emax"].real(emax)
    sim.run()

    # Generate an energy binning
    e_min   = gammalib.GEnergy(emin, "TeV")
    e_max   = gammalib.GEnergy(emax, "TeV")
    ebounds = gammalib.GEbounds(10, e_min, e_max)

    # Generate spectral points
    spectrum = obsutils.spectrum(sim.obs(), "Crab", ebounds)
	
    # Return spectrum
    return spectrum


# ============= #
# Plot spectrum #
# ============= #
def plot_spectrum(spectrum):
    """
    Plot spectrum.
    """
    # Create figure
    plt.figure(1)
    plt.title("Crab spectrum")

    # Plot spectrum
    plt.loglog(spectrum['energy']['value'], \
               spectrum['flux']['value'], 'ro', label='Crab')
    plt.errorbar(spectrum['energy']['value'], \
                 spectrum['flux']['value'], \
                 spectrum['flux']['ed_value'], ecolor='r')

    # Put labels
    plt.xlabel("Energy ("+spectrum['energy']['unit']+")")
    plt.ylabel("Flux ("+spectrum['flux']['unit']+")")
    plt.legend(loc="lower left")

    # Show spectrum
    plt.show()

    # Return
    return
    

#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    Generate a Crab spectrum.
    """
    # Generate spectrum
    spectrum = make_spectrum()

    # Plot spectrum
    if has_matplotlib:
        plot_spectrum(spectrum)
    