#! /usr/bin/env python
# ==========================================================================
# Shows histogram of residual map
#
# Copyright (C) 2018 Juergen Knoedlseder
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


# ============ #
# Create array #
# ============ #
def create_array(map,step_nx=10,step_ny=10,step_neng=1):
    """
    """
    # Initialise array of values
    values = []

    # Loop over sky map
    for ix in range(0,map.nx(),step_nx):
        for iy in range(0,map.ny(),step_ny):
            for ieng in range(0,map.nmaps(),step_neng):

                # Initialise value
                value = 0.0

                # Loop over sub block
                for ix_sub in range(step_nx):
                    ix_abs = ix + ix_sub
                    if ix_abs < map.nx():
                        for iy_sub in range(step_ny):
                            iy_abs = iy + iy_sub
                            if iy_abs < map.ny():
                                inx = ix_abs + iy_abs*map.nx()
                                for ieng_sub in range(step_neng):
                                    ieng_abs = ieng + ieng_sub
                                    if ieng_abs < map.nmaps():
                                        value += map(inx,ieng_abs)

                # If value differs from zero than append it
                if value != 0.0:
                    values.append(value)

    # Return values
    return values


# ================= #
# Plot residual map #
# ================= #
def plot_resmap(filename, nbins, plotfile):
    """
    Plot pull histogram

    Parameters
    ----------
    filename : str
        Residual map filename
    nbins : int
        Number of hisogram bins
    plotfile : str
        Plot filename
    """
    # Open residual map
    map = gammalib.GSkyMap(filename)

    # Create array of values
    values = create_array(map)

    # Create histogram
    _, bins, _ = plt.hist(values, nbins, range=[-4.0,4.0],
                          normed=True, facecolor='green')

    # Create expected distribution
    y = mlab.normpdf(bins, 0.0, 1.0)
    plt.plot(bins, y, 'r-', linewidth=2)

    # Set plot
    plt.xlabel('Pull')
    plt.ylabel('Arbitrary units')
    #plt.title(parname)
    plt.grid(True)

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ================= #
# Show residual map #
# ================= #
def show_resmap():
    """
    Show residual map
    """
    # Set usage string
    usage = 'show_resmap.py [-n bins] [-p plotfile] file'

    # Set default options
    options = [{'option': '-n', 'value': '50'},
               {'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    nbins    = int(options[0]['value'])
    plotfile = options[1]['value']

    # Plot residual map
    plot_resmap(args[0], nbins, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show residual map
    show_resmap()
