#! /usr/bin/env python
# ==========================================================================
# Makes spectrum and optionally plot it
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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
import ctools
import cscripts
try:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()


# =============== #
# Create spectrum #
# =============== #
def create_spectrum(datadir):
    """
    Simulate events and generate a source spectrum

    Returns
    -------
    datadir : str
        Data directory
    """
    # Set script parameters
    model_name  = datadir+'/crab.xml'
    caldb       = 'prod2'
    irf         = 'South_0.5h'
    ra          =   83.63
    dec         =   22.01
    rad_sim     =    3.0
    tstart      =    0.0
    tstop       = 1800.0
    emin        =    0.1
    emax        =  100.0
    enumbins    =     10

    # Simulate events
    sim = ctools.ctobssim()
    sim['inmodel'] = model_name
    sim['caldb']   = caldb
    sim['irf']     = irf
    sim['ra']      = ra
    sim['dec']     = dec
    sim['rad']     = rad_sim
    sim['tmin']    = tstart
    sim['tmax']    = tstop
    sim['emin']    = emin
    sim['emax']    = emax
    sim.run()

    # Create spectrum
    spec = cscripts.csspec(sim.obs())
    spec['srcname']  = 'Crab'
    spec['outfile']  = 'example_spectrum.fits'
    spec['expcube']  = 'NONE'
    spec['psfcube']  = 'NONE'
    spec['bkgcube']  = 'NONE'
    spec['edisp']    = False
    spec['emin']     = emin
    spec['emax']     = emax
    spec['enumbins'] = enumbins
    spec['ebinalg']  = 'LOG'
    spec.run()
    spec.save()

    # Get copy of spectrum
    spectrum = spec.spectrum().copy()

    # Return
    return spectrum


# ============= #
# Plot spectrum #
# ============= #
def plot_spectrum(spectrum, plotfile):
    """
    Plot spectrum

    Parameters
    ----------
    spectrum : `~gammalib.GFits`
        Spectrum FITS file
    plotfile : str
        Plot filename
    """
    # Extract columns from spectrum file
    table    = spectrum.table(1)
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

    # Get number of spectral points
    nrows = table.nrows()

    # Loop over rows of the file
    for row in range(nrows):

        # Get TS value, flux and flux error
        ts    = c_ts.real(row)
        flx   = c_flux.real(row)
        e_flx = c_eflux.real(row)

        # If the TS > 9 and the flux error is smaller than the flux then
        # take the fitted flux and flux error and append them to the
        # arrays
        if ts > 9.0 and e_flx < flx:

            # Add information
            energies.append(c_energy.real(row))
            flux.append(c_flux.real(row))
            ed_engs.append(c_ed.real(row))
            eu_engs.append(c_eu.real(row))
            e_flux.append(c_eflux.real(row))

        # ... otherwise take the upper limit and append it to the upper
        # limit array
        else:

            # Add information
            ul_energies.append(c_energy.real(row))
            ul_flux.append(c_upper.real(row))
            ul_ed_engs.append(c_ed.real(row))
            ul_eu_engs.append(c_eu.real(row))

    # Create figure
    plt.figure()
    plt.title('Crab spectrum')

    # Plot the spectrum 
    plt.loglog()
    plt.grid()

    # Plot error bars for flux points and upper limits for upper limit
    # points
    plt.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs],
                 fmt='ro')
    plt.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs],
                 yerr=1.0e-11, uplims=True, fmt='ro')

    # Add plot labels
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'E$^{2}$dN/dE (erg cm$^{-2}$ s$^{-1}$)')

    # Show spectrum or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ============= #
# Make spectrum #
# ============= #
def make_spectrum():
    """
    Make spectrum
    """
    # Set usage string
    usage = 'make_spectrum.py [-d datadir] [-p plotfile]'

    # Set default options
    options = [{'option': '-d', 'value': 'data'},
               {'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    datadir  = options[0]['value']
    plotfile = options[1]['value']

    # Generate spectrum
    spectrum = create_spectrum(datadir)

    # Plot spectrum
    plot_spectrum(spectrum, plotfile)

    # Return
    return

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Make spectrum
    make_spectrum()
