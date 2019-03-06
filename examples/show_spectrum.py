#! /usr/bin/env python
# ==========================================================================
# Display spectrum generated by csspec
#
# Copyright (C) 2015-2019 Juergen Knoedlseder
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
try:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()
import gammalib
import cscripts


# ============= #
# Plot spectrum #
# ============= #
def plot_spectrum(filename, plotfile):
    """
    Plot spectrum

    Parameters
    ----------
    filename : str
        Name of spectrum FITS file
    plotfile : str
        Plot file name
    """
    # Read spectrum file    
    fits     = gammalib.GFits(filename)
    table    = fits.table(1)
    c_energy = table['Energy']
    c_ed     = table['ed_Energy']
    c_eu     = table['eu_Energy']
    c_flux   = table['Flux']
    c_eflux  = table['e_Flux']
    c_ts     = table['TS']
    c_upper  = table['UpperLimit']

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

        # Get Test Statistic, flux and flux error
        ts    = c_ts.real(row)
        flx   = c_flux.real(row)
        e_flx = c_eflux.real(row)

        # If Test Statistic is larger than 9 and flux error is smaller than
        # flux then append flux plots ...
        if ts > 9.0 and e_flx < flx:
            energies.append(c_energy.real(row))
            flux.append(c_flux.real(row))
            ed_engs.append(c_ed.real(row))
            eu_engs.append(c_eu.real(row))
            e_flux.append(c_eflux.real(row))

        # ... otherwise append upper limit
        else:
            ul_energies.append(c_energy.real(row))
            ul_flux.append(c_upper.real(row))
            ul_ed_engs.append(c_ed.real(row))
            ul_eu_engs.append(c_eu.real(row))

    # Set upper limit errors
    yerr = [0.6 * x for x in ul_flux]

    # Plot the spectrum 
    plt.figure()
    plt.loglog()
    plt.grid()
    plt.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs],
                 fmt='ro')
    plt.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs],
                 yerr=yerr, uplims=True, fmt='ro')
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)')

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ============= #
# Show spectrum #
# ============= #
def show_spectrum():
    """
    Show spectrum
    """
    # Set usage string
    usage = 'show_spectrum.py [-p plotfile] [file]'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Show spectrum
    plot_spectrum(args[0], plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show spectrum
    show_spectrum()
