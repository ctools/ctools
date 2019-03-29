#! /usr/bin/env python
# ==========================================================================
# Shows source variability evolution created with ctfindvar
#
# Copyright (C) 2018 Simon Bonnefoy
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


# ========================== #
# Plot variability evolution #
# ========================== #
def plot_variability_evolution(filename, plotfile, srcname):
    """
    Plot variability evolution 

    Parameters
    ----------
    filename : str
        Name of spectrum FITS file
    plotfile : str
        Plot file name
    srcname : str
        Source name
    """
    # Read variability  file    
    fits    = gammalib.GFits(filename)
    table   = fits.table(1)
    c_gti_start  = table['TSTART']
    c_gti_stop   = table['TSTOP']
    c_sig        = table[srcname]

    # Convert the c_sig object to list
    table_c_sig = []
    table_c_gti = []
    for k in c_sig:
        table_c_sig.append(k)
    for j, k in zip(c_gti_start, c_gti_stop):
        table_c_gti.append(0.5*(j+k))

    # Initializing the plots
    axarr = []
    f     = plt.figure(figsize=(12,7))
    ax1   = f.add_subplot(121)
    ax2   = f.add_subplot(122)
    axarr.append(ax1)
    axarr.append(ax2)

    # Create first plot
    axarr[0].plot(table_c_gti, table_c_sig)
    axarr[0].set_title('Significance evolution')
    axarr[0].set_xlabel('MJD')
    axarr[0].set_ylabel('Significance')

    # Create second plot
    axarr[1].hist(table_c_sig, 20, range=[min(table_c_sig), max(table_c_sig)],
                  histtype='stepfilled', facecolor='g', alpha=0.75)
    axarr[1].set_title('Significance distribution')
    axarr[1].set_xlabel('Significance')
    axarr[1].set_ylabel('Counts')
    
    # Optionally save result in file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ========================== #
# Show variability evolution #
# ========================== #
def show_variability_evolution():
    """
    Show variability 
    """
    # Set usage string
    usage = 'show_variability_evolution.py [-p plotfile] [-s source_name] [file]'

    # Set default options
    options = [{'option': '-p', 'value': ''},
               {'option': '-s', 'value': 'MAXSIGPIXEL'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']
    srcname  = options[1]['value']

    # Show variability evolution
    plot_variability_evolution(args[0], plotfile, srcname)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show variability evolution 
    show_variability_evolution()
