#! /usr/bin/env python
# ==========================================================================
# Shows the distribution of significances in a given significance map.
#
# Copyright (C) 2011-2017 Andreas Specovius
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
    import matplotlib.mlab as mlab
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

try:
    import scipy.optimize as optimize
except ImportError:
    print('This script needs the "scipy" module')
    sys.exit()


#
#
#
# What do I need for plotting:
#
#- significance map itself [mappath]
# Provided as a file path guiding to a GSkyMap in FITS format
# No default
#
#- included regions (for shrinking data space i.e. for cutting certain edges or at a certain radius) [includedreg]
# Provided as a file path to a FITS WCS map or DS9 regions file
# Default should be None, which means all is allowed
#  
#- excluded regions [excludedreg]
# Provided as a file path to a FITS WCS map or DS9 regions file
# Default should be None, which means there are no excluded regions. In this case no bkg will be plotted and fitted
#
#- binning parameters for the histogram: [nbins], [sign_min], [sign_max]
# Provided as floats
# Default should be -5 to 8 with 130 bins (131 edges)
#
#- title of the histogram [title]
# Provided as string
# Default should be an empty string
#
#- Path to plotted file [plotfile]
#
#


# ================ #
# Helper functions #
# ================ #
def gaussian(x, norm, xmean, sigma):
    """
    Compute function values for gaussian function.
    """
    return norm*np.exp(-(x-xmean)**2/float(2*sigma**2))


def skymap_to_numpy_ndarray(skymap):
    """
    Return numpy.ndarray of pixel values of a skymap.
    """
    # Get map size
    npix  = skymap.npix()

    # Initialise out array
    array = np.empty(npix)

    # Loop over map pixels
    for ipix in range(npix):
        array[ipix] = skymap[ipix]

    # Return
    return array


def read_regions(fpath, init_map):
    """
    Read a regions WCS FITS or ds9 file and return a mask WCS map.

    Parameters
    ----------
    fpath : str
        Path to regions file
    init_map : GSkyMap
        Sky map representing fov of interest. Used for initialisation of
        regions map.
    """

    # Create starter map for user fov filled with zeros
    regions_map = init_map.copy()
    regions_map *= 0

    # Take care about ds9 region files
    if fpath.lower().endswith('.reg'):

        # Read ds9 regions
        ds9_regions = gammalib.GSkyRegions(fpath)

        # Loop over regions
        for reg in ds9_regions:

            # Make map from region
            ds9_map = gammalib.GSkyRegionMap(reg)

            # Add ds9 region map to global map
            regions_map += ds9_map.map()

    # Take care about fits wcs files
    elif fpath.lower().endswith('.fits'):

        # Read wcs regions map
        wcs_regions = gammalib.GSkyMap(fpath)

        # Add wcs regions map to global map
        regions_map += wcs_regions

    else:
        raise RuntimeError("Invalid regions file detected. Please provide " +
                           "a valid ds9 or FITS WCS regions file.")

    # Return
    return regions_map


# ============================== #
# Plot significance distribution #
# ============================== #
def plot_significance_distribution(mappath, nbins, sigma_min, sigma_max,
                                   includedreg, excludedreg, title):
    """
    Plot the significance distribution and return instance of pyplot
    figure and axis.

    Parameters
    ----------
    mappath : str
        Path to significance map
    nbins : int
        Number of bins in histogram
    sigma_min : float
        Lower limit of the x axis in the plot
    sigma_max : float
        Upper limit of the x axis in the plot
    includedreg : str
        Path to global included region file
    excludedreg : str
        Path to global excluded region file
    title : str
        Title of the plot

    Figure and axis objects are returned.
    """

    # Read significance map
    # =====================

    # Read map
    significance_map = gammalib.GSkyMap(mappath+'[SIGNIFICANCE]')


    # Read included regions
    # =====================

    # Set availability to false
    inclusion_available = False

    # Read in case provided
    if len(includedreg) > 0:
        # Read regions map
        regions_map = read_regions(includedreg, init_map=significance_map)

        # Make numpy boolean mask
        mask_include = skymap_to_numpy_ndarray(regions_map).astype(bool)

        # Enable masked map
        inclusion_available = True


    # Read excluded regions
    # =====================

    # Set availability to false
    src_exclusion_available = False

    # Read in case provided
    if len(excludedreg) > 0:
        # Read regions file
        regions_map = read_regions(excludedreg, init_map=significance_map)

        # Make numpy boolean mask
        mask_exclude = skymap_to_numpy_ndarray(regions_map).astype(bool)

        # Enable masked map
        src_exclusion_available = True


    # Convet maps to flat arrays and apply spatial cuts
    # =================================================

    # significance data to flat array
    data = skymap_to_numpy_ndarray(significance_map)

    # Initialise dummy mask allowing all pixels
    mask = np.ones(significance_map.npix()).astype(bool)

    # Apply inclusion region spatial cuts
    if inclusion_available:
        mask &= mask_include
        data = data[mask]

    # Apply exclusion region spatial cuts
    if src_exclusion_available:
        data_without_sources = skymap_to_numpy_ndarray(significance_map)
        mask &= ~mask_exclude
        data_without_sources = data_without_sources[mask]


    # Compute bin edges to user parameters
    # ====================================

    # Compute bin edges, hence use nbins+1
    bin_edges = np.linspace(sigma_min, sigma_max, nbins+1)


    # Begin with plotting
    # ===================

    # Create figure
    fig = plt.figure()
    ax  = fig.gca()

    # Draw significance histogram for full fov
    y, _, _ = ax.hist(data, bins=bin_edges, histtype='step', color='k', label="significance")

    # Draw the default normal distribution. Scale to fit maximum of histogram
    ax.plot(bin_edges, gaussian(bin_edges, y.max(),0,1), 'b-')
    msg = 'mean: $0$\nwidth: $1$'
    ax.text(0.98, 0.80, msg, ha='right', va='top', bbox=dict(edgecolor='blue', facecolor='white'), transform=ax.transAxes)

    # If exclusion was provided also handle significanes with excluded src data
    if src_exclusion_available:

        # Draw histogram of significances for data without excluded sources
        y, _, _ = ax.hist(data_without_sources, bins=bin_edges, histtype='step',
                          color='k', alpha=0.5, label="significance without exclusions")

        # Fit a normal distribution to the binned masked significance data
        pars        = [1.0, 0.0, 1.0]
        bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

        pars_opt, pars_cov = optimize.curve_fit(gaussian, xdata=bin_centers,
                                                ydata=y, p0=pars)

        # Compute parameters fit errors
        pars_err = np.sqrt(np.diag(pars_cov))

        # Log fit result
        print('Fit of normal distribution to excluded source significance data ' +\
              'results in:  norm=%f+=%f  xmean=%f+-%f sigma=%f+-%f' %
              (pars_opt[0], pars_err[0], pars_opt[1], pars_err[1], pars_opt[2], pars_err[2])
             )

        # Draw text box
        msg = 'mean: $%.3f\pm%.3f$\nwidth: $%.3f\pm%.3f$' % (pars_opt[1], pars_err[1], pars_opt[2], pars_err[2])
        ax.text(0.98, 0.95, msg, ha='right', va='top', bbox=dict(edgecolor='red', facecolor='white'), transform=ax.transAxes)

        # Plot the normal (2)
        ax.plot(bin_edges, gaussian(bin_edges, *(pars_opt)), 'r-')

    # Configure the plot
    ax.set_ylim(0.5, ax.get_ylim()[1]*2.0)
    ax.set_xlim(sigma_min, sigma_max)
    ax.set_title(title)
    ax.set_xlabel("Significance")
    ax.set_ylabel("Entries")
    ax.grid()
    ax.set_yscale('log')

    # Return
    return fig, ax


# ============================== #
# Show significance distribution #
# ============================== #
def show_significance_distribution():
    """
    Plot and show significance distribution
    """
    # Set usage string
    usage = 'show_significance_distribution.py [--nbins bins] ' +\
            '[--smin sigma_min] [--smax sigma_max] [--include regfile] ' +\
            '[--exclude regfile] [--title title] [-p plotfile] mapfile'

    # Set default options
    options = [{'option': '--nbins',    'value': '130'},
               {'option': '--smin',     'value':  '-5'},
               {'option': '--smax',     'value':   '8'},
               {'option': '--include',  'value':   ''},
               {'option': '--exclude',  'value':   ''},
               {'option': '--title',    'value':   ''},
               {'option': '-p',         'value':   ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Print help message
    if args[0]=='--help' or args[0]=='-h':
        print(usage)
        sys.exit()

    # Extract script parameters from options
    nbins       = int(options[0]['value'])
    sigma_min   = int(options[1]['value'])
    sigma_max   = int(options[2]['value'])
    includedreg = options[3]['value']
    excludedreg = options[4]['value']
    title       = options[5]['value']
    plotfile    = options[6]['value']

    # Plot significance distribution histogram
    fig, ax = plot_significance_distribution(args[0], nbins, sigma_min, sigma_max,
                                             includedreg, excludedreg, title)

    # Show or save the plot
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Plot and show significance distribution histogram
    show_significance_distribution()
