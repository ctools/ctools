#! /usr/bin/env python
# ==========================================================================
# Display spectrum generated by csspec
#
# Copyright (C) 2015-2020 Juergen Knoedlseder
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
import gammalib
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


# ===================================== #
# Extract the spectrum info from a file #
# ===================================== #
def get_spectrum_file(filename):
    """
    Extract the spectrum info from a file for plotting

    Parameters
    ----------
    filename : str
        Name of spectrum FITS file

    Returns
    -------
    spec : dict
        Python dictionary defining spectral plot parameters
    """
    # Read spectrum file    
    fits = gammalib.GFits(filename)

    # Return dictionary
    return get_spectrum_fits(fits)


# ============================================= #
# Extract the spectrum info from a GFits object #
# ============================================= #
def get_spectrum_fits(fits):
    """
    Extract the spectrum info from a GFits object

    Parameters
    ----------
    fits : `~gammalib.GFits`
        Spectral GFits object
    
    Returns
    -------
    spec : dict
        Python dictionary defining spectral plot parameters
    """
    # Read spectrum objects
    table    = fits.table(1)
    c_energy = table['e_ref']
    c_ed     = table['e_min']
    c_eu     = table['e_max']
    c_flux   = table['ref_e2dnde']
    c_norm   = table['norm']
    c_eflux  = table['norm_err']
    c_upper  = table['norm_ul']
    c_ts     = table['ts']

    # Initialise arrays to be filled
    spec = {
        'energies'    : [],
        'flux'        : [],
        'ed_engs'     : [],
        'eu_engs'     : [],
        'e_flux'      : [],
        'ul_energies' : [],
        'ul_ed_engs'  : [],
        'ul_eu_engs'  : [],
        'ul_flux'     : [],
        'yerr'        : [],
        'loglike'     : [],
        'engs_scan'   : [],
        'e2dnde_scan' : [],
        'dll_scan'    : []
    }

    # Determine if we can load the delta-log-likelihood profiles
    has_sedtype = table.has_card('SED_TYPE')
    load_dll    = False
    if has_sedtype:
        seds = table.card('SED_TYPE').string().split(',')

        # Load likelhood columns if present
        if seds[0] == 'likelihood':
            load_dll    = True
            c_loglike   = table['loglike']
            c_norm_scan = table['norm_scan']
            c_dll_scan  = table['dloglike_scan']

    # Loop over rows of the file
    nrows = table.nrows()
    for row in range(nrows):

        # Get Test Statistic, flux and flux error
        ts    = c_ts.real(row)
        norm  = c_norm.real(row)
        flx   = norm * c_flux.real(row)
        e_flx = flx * c_eflux.real(row)

        # If Test Statistic is larger than 9 and flux error is smaller than
        # flux then append flux plots ...
        if ts > 9.0 and e_flx < flx:
            spec['energies'].append(c_energy.real(row))
            spec['flux'].append(flx)
            spec['ed_engs'].append(c_ed.real(row))
            spec['eu_engs'].append(c_eu.real(row))
            spec['e_flux'].append(e_flx)

        # ... otherwise append upper limit
        else:
            spec['ul_energies'].append(c_energy.real(row))
            spec['ul_flux'].append(flx*c_upper.real(row))
            spec['ul_ed_engs'].append(c_ed.real(row))
            spec['ul_eu_engs'].append(c_eu.real(row))

    # Load the delta log-likelihood values
    if load_dll:
        spec['e2dnde_scan'] = np.zeros((nrows, c_norm_scan.elements(0)))
        spec['dll_scan']    = np.zeros((nrows, c_norm_scan.elements(0)))

        for row in range(nrows):
            loglike = c_loglike.real(row)
            spec['loglike'].append(loglike)
            spec['engs_scan'].append(c_energy.real(row) - c_ed.real(row))
            for col in range(c_norm_scan.elements(row)):
                spec['e2dnde_scan'][row,col] = c_norm_scan.real(row,col) * c_flux.real(row)
                spec['dll_scan'][row,col] = c_dll_scan.real(row,col)
        spec['engs_scan'].append(c_energy.real(nrows-1)+c_eu.real(nrows-1))

    # Set upper limit errors
    spec['yerr'] = [0.6 * x for x in spec['ul_flux']]

    # Return dictionary
    return spec


# ========================= #
# Plot delta log-likelihood #
# ========================= #
def plot_dloglike(spec):
    """
    Plot delta log-likelihood

    Parameters
    ----------
    spec : dict
        Python dictionary defining spectral plot parameters
    """
    # Establish the bounds on the x,y plot
    ymin,ymax = plt.gca().get_ylim()
    ymin      = np.log10(ymin)
    ymax      = np.log10(ymax)
    steps     = 1000
    ebins     = len(spec['e2dnde_scan'])
    y_bounds  = np.linspace(ymin, ymax, steps+1)

    # Compute the logarithmic center of each bin
    fluxpnts = 10 ** ((y_bounds[1:]+y_bounds[:-1]) / 2.0)
    
    # Scale the boundaries on the grid
    y_bounds = 10 ** y_bounds

    # Interpolate the dlogLike values
    dll_hist = [[]]*ebins
    for ebin in range(ebins):
        dll_hist[ebin] = np.interp(fluxpnts, spec['e2dnde_scan'][ebin], 
                                             spec['dll_scan'][ebin])

    # Transpose the dll_hist array so x=energy, y=dnde
    dll_hist = [[row[i] for row in dll_hist] for i in range(steps)]

    # Plot the likelihood profiles
    plt.pcolormesh(np.array(spec['engs_scan']), y_bounds,
                   np.array(dll_hist), vmin=-10.0, vmax=0.0,
                   cmap='Reds', alpha=0.9)
    cbar = plt.colorbar()
    cbar.set_label(r'$\Delta$ log-likelihood')
    plt.grid()

    # Return
    return


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
    # Get spectrum parameters
    spec = get_spectrum_file(filename)

    # Create the plot
    plt.figure()
    plt.loglog()
    plt.grid()
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)')

    # Plot the spectrum
    plt.errorbar(spec['energies'], spec['flux'], 
                 yerr=spec['e_flux'], xerr=[spec['ed_engs'], spec['eu_engs']],
                 fmt='ro')

    # Plot upper limits
    if len(spec['ul_energies'])	> 0:
        plt.errorbar(spec['ul_energies'], spec['ul_flux'], yerr=spec['yerr'],
                     xerr=[spec['ul_ed_engs'], spec['ul_eu_engs']],
                     uplims=True, fmt='ro')

    # Plot log-likelihood profiles
    if len(spec['dll_scan']) > 0:
        plot_dloglike(spec)

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
