#! /usr/bin/env python
# ==========================================================================
# Shows the photon spectrum of a model
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
import gammalib
import cscripts
try:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()


# ======================= #
# Plot spectral component #
# ======================= #    
def plot_spectrum(model, plotfile, emin=0.01, emax=100.0, enumbins=100):
    """
    Plot spectral model component

    Parameters
    ----------
    model : `~gammalib.GModel`
        Model
    plotfile : str
        Name of plot file
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    enumbins : integer, optional
        Number of energy bins
    """
    # Get spectral component
    spectrum = model.spectral()

    # Setup energy axis
    e_min   = gammalib.GEnergy(emin, 'TeV')
    e_max   = gammalib.GEnergy(emax, 'TeV')
    ebounds = gammalib.GEbounds(enumbins, e_min, e_max)

    # Setup lists of x and y values
    x   = []
    y   = []
    for i in range(enumbins):
        energy = ebounds.elogmean(i)
        value  = spectrum.eval(energy)
        x.append(energy.TeV())
        y.append(value)

    # Show spectrum
    plt.figure()
    plt.title(model.name()+' ('+spectrum.type()+')')
    plt.loglog()
    plt.grid()
    plt.loglog(x, y, color='red')
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'dN/dE (ph s$^{-1}$ cm$^{-2}$ MeV$^{-1}$)')

    # Show spectrum or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# =================== #
# Show model spectrum #
# =================== #
def show_model():
    """
    Show model spectrum
    """
    # Set usage string
    usage = 'show_model.py [-n name] [-p plotfile] file'

    # Set default options
    options = [{'option': '-n', 'value': ''},
               {'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    name     = options[0]['value']
    plotfile = options[1]['value']
    if len(name) == 0:
        name = 0

    # Read models XML file
    models = gammalib.GModels(args[0])

    # Extract relevant model
    model = models[name]

    # Plot spectrum
    plot_spectrum(model, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show model spectrum
    show_model()
