#! /usr/bin/env python
# ==========================================================================
# Shows butterfly diagram created with ctbutterfly
#
# Copyright (C) 2014-2022 Michael Mayer
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


# =================== #
# Read butterfly data #
# =================== #
def read_butterfly(filename):
    """
    Read butterfly data

    Parameters
    ----------
    filename : str
        Name of FITS or CSV file

    Returns
    -------
    butterfly : dict
        Python dictionary defining butterfly plot and best fit spectrum
    """
    # Read butterfly data dependent on file type
    fname = gammalib.GFilename(filename)
    if fname.is_fits():
        butterfly = read_butterfly_fits(filename)
    else:
        butterfly = read_butterfly_csv(filename)

    # Return butterfly dictionary
    return butterfly


# ================================= #
# Read butterfly data from CSV file #
# ================================= #
def read_butterfly_csv(filename):
    """
    Read butterfly data from CSV file

    Parameters
    ----------
    filename : str
        Name of CSV file

    Returns
    -------
    butterfly : dict
        Python dictionary defining butterfly plot and best fit spectrum
    """
    # Initialise arrays to be filled
    butterfly = {'butterfly_x' : [],
                 'butterfly_y' : [],
                 'line_x'      : [],
                 'line_y'      : []}

    # Open butterfly file
    csv = gammalib.GCsv(filename)

    # Loop over rows of the file
    nrows = csv.nrows()
    for row in range(nrows):

        # Get conversion coefficient
        conv = csv.real(row,0) * csv.real(row,0) * gammalib.MeV2erg

        # Compute upper edge of confidence band
        butterfly['butterfly_x'].append(csv.real(row,0)/1.0e6) # TeV
        butterfly['butterfly_y'].append(csv.real(row,2)*conv)

        # Set line values
        butterfly['line_x'].append(csv.real(row,0)/1.0e6) # TeV
        butterfly['line_y'].append(csv.real(row,1)*conv)

    # Loop over the rows backwards to compute the lower edge of the
    # confidence band
    for row in range(nrows):
        index     = nrows - 1 - row
        conv      = csv.real(index,0) * csv.real(index,0) * gammalib.MeV2erg
        low_error = max(csv.real(index,3)*conv, 1e-26)
        butterfly['butterfly_x'].append(csv.real(index,0)/1.0e6)
        butterfly['butterfly_y'].append(low_error)

    # Return butterfly dictionary
    return butterfly


# ================================== #
# Read butterfly data from FITS file #
# ================================== #
def read_butterfly_fits(filename):
    """
    Read butterfly data from FITS file

    Parameters
    ----------
    filename : str
        Name of FITS file

    Returns
    -------
    butterfly : dict
        Python dictionary defining butterfly plot and best fit spectrum
    """
    # Initialise arrays to be filled
    butterfly = {'butterfly_x' : [],
                 'butterfly_y' : [],
                 'line_x'      : [],
                 'line_y'      : []}

    # Open FITS file
    fits = gammalib.GFits(filename)

    # Get sensitivity table
    table = fits.table('BUTTERFLY')

    # Get relevant columns
    c_energy        = table['ENERGY']
    c_intensity     = table['INTENSITY']
    c_intensity_min = table['INTENSITY_MIN']
    c_intensity_max = table['INTENSITY_MAX']

    # Fill vectors
    nrows = table.nrows()
    for row in range(nrows):

        # Exclude zero intensities
        if c_intensity[row] <= 0.0:
            continue

        # Get conversion coefficient TeV -> erg
        conv = c_energy[row] * c_energy[row] * 1.0e6 * gammalib.MeV2erg

        # Compute upper edge of confidence band
        butterfly['butterfly_x'].append(c_energy[row])
        butterfly['butterfly_y'].append(c_intensity_max[row] * conv)

        # Set line values
        butterfly['line_x'].append(c_energy[row])
        butterfly['line_y'].append(c_intensity[row] * conv)

    # Loop over the rows backwards to compute the lower edge of the
    # confidence band
    for row in range(nrows-1,-1,-1):
        if c_intensity[row] <= 0.0:
            continue
        conv      = c_energy[row] * c_energy[row] * 1.0e6 * gammalib.MeV2erg
        low_error = max(c_intensity_min[row] * conv, 1e-26)
        butterfly['butterfly_x'].append(c_energy[row])
        butterfly['butterfly_y'].append(low_error)

    # Return butterfly dictionary
    return butterfly


# ============== #
# Plot butterfly #
# ============== #
def plot_butterfly(filename, plotfile):
    """
    Plot butterfly diagram

    Parameters
    ----------
    filename : str
        Butterfly CSV file
    plotfile : str
        Plot filename
    """
    # Read butterfly data
    butterfly = read_butterfly(filename)

    # Plot the butterfly and spectral line       
    plt.figure()
    plt.loglog()
    plt.grid()
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)')

    # Plot the butterfly and spectral line
    plt.plot(butterfly['line_x'],      butterfly['line_y'],
             color='black', ls='-')
    plt.fill(butterfly['butterfly_x'], butterfly['butterfly_y'],
             color='green', alpha=0.5)

    # Show spectrum or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ====================== #
# Show butterfly diagram #
# ====================== #
def show_butterfly():
    """
    Show butterfly diagram
    """
    # Set usage string
    usage = 'show_butterfly.py [-p plotfile] file'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Plot butterfly diagram
    plot_butterfly(filename=args[0], plotfile=plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show butterfly diagram
    show_butterfly()
