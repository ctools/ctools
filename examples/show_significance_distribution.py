#! /usr/bin/env python
# ==========================================================================
# Shows the distribution of significances in a given significance map.
#
# Copyright (C) 2018-2019 Andreas Specovius
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
# --------------------------------------------------------------------------
# What do I need for plotting:
#
# - significance map itself [mappath]
#   Provided as a file path guiding to a GSkyMap in FITS format
#   No default
#
# - included regions (for shrinking data space i.e. for cutting certain edges
#   or at a certain radius) [includedreg]
#   Provided as a file path to a FITS WCS map or DS9 regions file
#   Default should be None, which means all is allowed
#  
# - excluded regions [excludedreg]
#   Provided as a file path to a FITS WCS map or DS9 regions file
#   Default should be None, which means there are no excluded regions. In this
#   case no bkg will be plotted and fitted
#
# - binning parameters for the histogram: [nbins], [sign_min], [sign_max]
#   Provided as floats
#   Default should be -5 to 8 with 130 bins (131 edges)
#
# - title of the histogram [title]
#   Provided as string
#   Default should be an empty string
#
# - Path to plotted file [plotfile]
# ==========================================================================
import sys
import math
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


# ======================= #
# Gaussian function class #
# ======================= #
class gaussian(gammalib.GPythonOptimizerFunction):

    # Constructor
    def __init__(self, x_vals, y_vals):
    
        # Call base class constructor
        gammalib.GPythonOptimizerFunction.__init__(self)
        
        # Set eval method
        self._set_eval(self.eval)

        # Set data
        self._x_vals = x_vals
        self._y_vals = y_vals

    # Methods
    def eval(self):
        """
        Evaluate function
        """
        # Recover parameters
        pars  = self._pars()
        norm  = pars[0].value()
        mean  = pars[1].value()
        sigma = pars[2].value()

        # Evaluate function values
        y = [norm * math.exp(-0.5*(x-mean)**2/(sigma**2)) for x in self._x_vals]

        # Compute weights (1/sqrt(y))
        weight = []
        for val in self._y_vals:
            if val > 0.0:
                weight.append(1.0/val)
            else:
                weight.append(0.0)
        #weight = [1.0/val if val > 0.0 else 0.0 for val in self._y_vals]

        # Compute Chi Square
        value = 0.0
        for i in range(len(self._x_vals)):
            arg    = self._y_vals[i] - y[i]
            value += arg * arg * weight[i]
        
        # Evaluate gradient and curvature
        sigma2 = sigma  * sigma
        sigma3 = sigma2 * sigma
        for i in range(len(self._x_vals)):

            # Evaluate function gradients
            dx     = self._x_vals[i] - mean
            dnorm  = y[i]         / norm   * pars[0].scale()
            dmean  = y[i] * dx    / sigma2 * pars[1].scale()
            dsigma = y[i] * dx**2 / sigma3 * pars[2].scale()

            # Setup gradient vector
            arg                 = (self._y_vals[i] - y[i]) * weight[i]
            self.gradient()[0] -= arg * dnorm
            self.gradient()[1] -= arg * dmean
            self.gradient()[2] -= arg * dsigma

            # Setup curvature matrix
            self.curvature()[0,0] +=  dnorm  * dnorm   * weight[i]
            self.curvature()[0,1] +=  dnorm  * dmean   * weight[i]
            self.curvature()[0,2] +=  dnorm  * dsigma  * weight[i]
            self.curvature()[1,0] +=  dmean  * dnorm   * weight[i]
            self.curvature()[1,1] +=  dmean  * dmean   * weight[i]
            self.curvature()[1,2] +=  dmean  * dsigma  * weight[i]
            self.curvature()[2,0] +=  dsigma * dnorm   * weight[i]
            self.curvature()[2,1] +=  dsigma * dmean   * weight[i]
            self.curvature()[2,2] +=  dsigma * dsigma  * weight[i]

        # Set value
        self._set_value(value)
        
        # Return
        return


# ================ #
# Helper functions #
# ================ #
def skymap_to_numpy_ndarray(skymap):
    """
    Return numpy.ndarray of pixel values of a skymap.

    Parameters
    ----------
    skymap : `~gammalib.GSkyMap`
        Sky map.

    Returns
    -------
    array : `~numpy.array`
        Numpy array
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


def regions_to_map(regions, init_map):
    """
    Generate a regions map based on an input skymap and a set of regions

    Parameters
    ----------
    region : `~gammalib.GSkyRegions`
        List of regions to be filled into a sky map
    init_map : `~gammalib.GSkyMap`
        Sky map representing fov of interest. Used for initialisation of
        regions map.

    Returns
    -------
    map : `~gammalib.GSkyMap`
        Sky map.
    """
    # Create starter map for user fov filled with zeros
    regions_map  = init_map.copy()
    regions_map *= 0.0

    # Loop over regions
    for reg in regions:

        # Make map from region
        ds9_map = gammalib.GSkyRegionMap(reg)

        # Add ds9 region map to global map
        regions_map += ds9_map.map()

    return regions_map


def read_regions(fpath):
    """
    Read a regions WCS FITS or ds9 file and return a GSkyRegions container.

    Parameters
    ----------
    fpath : str
        Path to regions file

    Returns
    -------
    regions : `~gammalib.GSkyRegions`
        Collection of regions specified in fpath
    """
    regions = gammalib.GSkyRegions()

    # Take care about ds9 region files
    if fpath.lower().endswith('.reg'):

        # Read ds9 regions
        regions = gammalib.GSkyRegions(fpath)

    # Take care about fits WCS files
    elif gammalib.GFilename(fpath).is_fits():

        # Read wcs regions map
        wcs_regions = gammalib.GSkyRegionMap(fpath)

        # Add wcs regions map to global regions
        regions.append(wcs_regions)

    else:
        raise RuntimeError('Invalid regions file detected. Please provide ' +
                           'a valid ds9 or FITS WCS regions file.')

    # Return
    return regions


def linspace(minval, maxval, nbins):
    """
    Return a linearly spaced array
    
    Parameters
    ----------
    minval : float
        Minimum value
    maxval : float
        Maximum value
    nbins : int
        Number of bins

    Returns
    -------
    bins : list of float
        Linearly spaced array
    """
    # Compute bin width
    binwidth = (maxval - minval) / float(nbins-1)

    # Compute bins
    bins = [minval + binwidth*i for i in range(nbins)]

    # Return bins
    return bins


# ============================== #
# Plot significance distribution #
# ============================== #
def plot_significance_distribution(mappath, nbins, sigma_min, sigma_max,
                                   includedreg, excludedreg, title,
                                   plotfile):
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
    plotfile : str
        Name of file for plotting
    """

    # Read significance map
    significance_map = gammalib.GSkyMap(mappath+'[SIGNIFICANCE]')

    # Read inclusion region file
    include_regs = None
    if len(includedreg) > 0:
        include_regs = read_regions(includedreg)

    # Read exclusion region file
    exclude_regs = None
    if len(excludedreg) > 0:
        exclude_regs = read_regions(excludedreg)

    # Generate the plot
    plot_significance(significance_map, nbins, sigma_min, sigma_max,
                      includedregs=include_regs, excludedregs=exclude_regs,
                      title=title, plotfile=plotfile)

    return


def plot_significance(sigmap, nbins, sigma_min, sigma_max,
                      includedregs=None, excludedregs=None, 
                      title='', plotfile=''):
    """
    Plot the significance distribution and return instance of pyplot
    figure and axis.

    Parameters
    ----------
    sigmap : `~gammalib.GSkyMap`
        Significance sky map
    nbins : int
        Number of bins in histogram
    sigma_min : float
        Lower limit of the x axis in the plot
    sigma_max : float
        Upper limit of the x axis in the plot
    includedregs : `~gammalib.GSkyRegions`
        Container of inclusion regions
    excludedregs : `~gammalib.GSkyRegions`
        Container for exclusion regions
    title : str
        Title of the plot
    plotfile : str
        Name of file for plotting
    """

    # Read included regions
    inclusion_available = False
    if includedregs is not None:
        regions_map         = regions_to_map(includedregs, init_map=sigmap)
        mask_include        = skymap_to_numpy_ndarray(regions_map).astype(bool)
        inclusion_available = True

    # Read excluded regions
    src_exclusion_available = False
    if excludedregs is not None:
        regions_map             = regions_to_map(excludedregs, init_map=sigmap)
        mask_exclude            = skymap_to_numpy_ndarray(regions_map).astype(bool)
        src_exclusion_available = True

    # Convert significance map to flat array
    data = skymap_to_numpy_ndarray(sigmap)

    # Initialise dummy mask allowing all pixels
    mask = np.ones(sigmap.npix()).astype(bool)

    # Apply inclusion region spatial cuts
    if inclusion_available:
        mask &= mask_include
        data = data[mask]

    # Apply exclusion region spatial cuts
    if src_exclusion_available:
        data_without_sources = skymap_to_numpy_ndarray(sigmap)
        mask &= ~mask_exclude
        data_without_sources = data_without_sources[mask]

    # Compute bin edges, hence use nbins+1
    bin_edges = linspace(sigma_min, sigma_max, nbins+1)

    # Create figure
    fig = plt.figure()
    ax  = fig.gca()

    # Draw significance histogram for full FOV
    y, _, _ = ax.hist(data, bins=bin_edges, histtype='step', color='k',
                      label='significance')

    # Draw the default normal distribution. Scale to fit maximum of histogram
    yvals = [y.max() * math.exp(-0.5*(x-0.0)**2/(1.0**2)) for x in bin_edges]
    ax.plot(bin_edges, yvals, 'b-')
    msg = 'mean: $0$\nwidth: $1$'
    ax.text(0.98, 0.80, msg, ha='right', va='top',
            bbox=dict(edgecolor='blue', facecolor='white'),
            transform=ax.transAxes)

    # If exclusion was provided then also plot significances with excluded
    # source data
    if src_exclusion_available:

        # Draw histogram of significances for data with excluded sources
        y, _, _ = ax.hist(data_without_sources, bins=bin_edges, histtype='step',
                          color='k', alpha=0.5,
                          label='significance without exclusions')

        # Set initial Gaussian parameters
        y_max = float(y.max())
        par1  = gammalib.GOptimizerPar('Norm',  y_max)
        par2  = gammalib.GOptimizerPar('Mean',  0.0)
        par3  = gammalib.GOptimizerPar('Sigma', 1.0)
        pars  = gammalib.GOptimizerPars()
        pars.append(par1)
        pars.append(par2)
        pars.append(par3)

        # Set fit function
        x   = [0.5*(bin_edges[i]+bin_edges[i+1]) for i in range(nbins)]
        fct = gaussian(x, y)

        # Optimize function and compute errors
        opt = gammalib.GOptimizerLM()
        opt.optimize(fct, pars)
        opt.errors(fct, pars)

        # Recover parameters and errors
        norm    = pars[0].value()
        e_norm  = pars[0].error()
        mean    = pars[1].value()
        e_mean  = pars[1].error()
        sigma   = pars[2].value()
        e_sigma = pars[2].error()

        # Log fit result
        print('Fit of normal distribution to excluded source significance data '+
              'results in:  norm=%f+=%f  xmean=%f+-%f sigma=%f+-%f' %
              (norm, e_norm, mean, e_mean, sigma, e_sigma))

        # Draw text box
        msg = 'mean: $%.3f\pm%.3f$\nwidth: $%.3f\pm%.3f$' % \
              (mean, e_mean, sigma, e_sigma)
        ax.text(0.98, 0.95, msg, ha='right', va='top',
                bbox=dict(edgecolor='red', facecolor='white'),
                transform=ax.transAxes)

        # Plot the normal
        yvals = [norm * math.exp(-0.5*(x-mean)**2/(sigma**2)) for x in bin_edges]
        ax.plot(bin_edges, yvals, 'r-')

    # Configure the plot
    ax.set_ylim(0.5, ax.get_ylim()[1]*2.0)
    ax.set_xlim(sigma_min, sigma_max)
    ax.set_title(title)
    ax.set_xlabel('Significance')
    ax.set_ylabel('Entries')
    ax.grid()
    ax.set_yscale('log')

    # Show or save the plot
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


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
    plot_significance_distribution(args[0], nbins, sigma_min, sigma_max,
                                   includedreg, excludedreg, title,
                                   plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Plot and show significance distribution histogram
    show_significance_distribution()
