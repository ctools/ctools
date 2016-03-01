#! /usr/bin/env python
# ==========================================================================
# This script shows a spectrum created with csspec
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

import gammalib
import sys
try:
    import matplotlib.pyplot as plt
except:
    sys.exit("This script needs matplotlib")


# ============ #
# Script entry #
# ============ #    
if __name__ == "__main__":
    """
    """

    # Check for given spectrum file
    if not len(sys.argv) == 2:
        sys.exit("Usage: show_spectrum.py spectrum.fits")

    # Read spectrum file    
    filename = sys.argv[1]
    fits     = gammalib.GFits(filename)
    table    = fits.table(1)
    c_energy = table["Energy"]
    c_ed     = table["ed_Energy"]
    c_eu     = table["eu_Energy"]
    c_flux   = table["Flux"]
    c_eflux  = table["e_Flux"]
    c_ts     = table["TS"]
    c_upper  = table["UpperLimit"]

    # Initialise arrays to be filled
    energies    = []
    flux        = []
    ed_engs     = []
    eu_engs     = []
    e_flux      = []
    ul_energies = []
    ul_ed_engs  = []
    ul_eu_engs  = []
    ul_flux     = []

    # Loop over rows of the file
    nrows = table.nrows()
    for row in range(nrows):

        # Get TS
        ts    = c_ts.real(row)
        flx   = c_flux.real(row)
        e_flx = c_eflux.real(row)

        # Switch
        if ts > 9.0 and e_flx < flx:

            # Add information
            energies.append(c_energy.real(row))
            flux.append(c_flux.real(row))
            ed_engs.append(c_ed.real(row))
            eu_engs.append(c_eu.real(row))
            e_flux.append(c_eflux.real(row))

        #
        else:

            # Add information
            ul_energies.append(c_energy.real(row))
            ul_flux.append(c_upper.real(row))
            ul_ed_engs.append(c_ed.real(row))
            ul_eu_engs.append(c_eu.real(row))

    # Plot the spectrum 
    plt.figure()
    plt.loglog()
    plt.grid()
    plt.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs], fmt='ro')
    plt.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs], yerr=1.0e-11, uplims=True, fmt='ro')
    plt.xlabel("Energy (TeV)")
    plt.ylabel(r"E dN/dE (erg cm$^{-2}$ s$^{-1}$)")    
    plt.show()
