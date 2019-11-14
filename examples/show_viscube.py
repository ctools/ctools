#! /usr/bin/env python
# ==========================================================================
# Show visibility cube
#
# Copyright (C) 2019 Juergen Knoedlseder
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
import matplotlib.pyplot as plt


# =========================== #
# Convert map into Aitoff map #
# =========================== #
def get_aitoff_map(filename):
    """
    Convert sky map file into an Aitoff map

    Parameters
    ----------
    filename : str
        File name

    Returns
    -------
    aitoff = `~gammalib.GSkyMap`
        Aitoff sky map
    """
    # Read sky map
    map = gammalib.GSkyMap(filename)

    # Integrate map over all zenith angles
    map.stack_maps()

    # Allocate Aitoff map
    aitoff = gammalib.GSkyMap('AIT', 'GAL', 0.0, 0.0, -1.0, 1.0, 324, 162)

    # Loop over all pixels in Aitoff map
    for i in range(aitoff.npix()):

        # Get sky direction
        try:
            # Get direction
            dir = aitoff.inx2dir(i)

            # Get value
            value = map(dir)

            # Set pixel value
            aitoff[i] = value

        except:
            continue

    # Return
    return aitoff


# ==================== #
# Plot visibility cube #
# ==================== #
def plot_viscube(ax, filename):
    """
    Plot visibility cube

    Parameters
    ----------
    ax : plt
        Subplot
    filename : str
        File name
    """
    # Get Aitoff map
    map = get_aitoff_map(filename)
    
    # Create array from skymap
    array = []
    for iy in range(map.ny()):
        row = []
        for ix in range(map.nx()):
            index = ix+iy*map.nx()
            value = map[index]
            row.append(value)
        array.append(row)

    # Show Aitoff projection
    c = ax.imshow(array,
                  extent=(gammalib.pi,-gammalib.pi,gammalib.pi/2,-gammalib.pi/2),
                  aspect=0.5)
    cbar = plt.colorbar(c, ax=ax, orientation='vertical', shrink=0.7)
    cbar.set_label('Visibility (hours)')
    ax.grid(True)
    ax.set_xlabel('Galactic longitude (deg)')
    ax.set_ylabel('Galactic latitude (deg)')
    ax.set_xticklabels([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])

    # Return
    return


# ======================= #
# Plot visibility results #
# ======================= #
def plot_visibility(ax, filename):
    """
    Plot visibility results

    Parameters
    ----------
    ax : plt
        Subplot
    filename : str
        File name
    """
    # Load results
    fits      = gammalib.GFits(filename)
    table     = fits.table('VISIBILITY')
    mjd       = table['MJD']
    dark_time = table['Darktime']

    # Create arrays
    nrows = table.nrows()
    time  = []
    dark  = []
    for row in range(nrows):
        time.append(mjd.real(row))
        dark.append(dark_time.real(row))

    # Plot arrays
    ax.plot(time, dark, 'r-')
    ax.set_xlabel('MJD (days)')
    ax.set_ylabel('Dark time (hours)')

    # Return
    return


# ==================== #
# Show visibility cube #
# ==================== #
def show_viscube():
    """
    Show visibility cube
    """
    # Set usage string
    usage = 'show_viscube.py [-t title] [-p plotfile] [filename]'

    # Set default options
    options = [{'option': '-p', 'value': ''},
               {'option': '-t', 'value': 'Visibility map'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']
    title    = options[1]['value']

    # Create figure and subplots
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(211, projection='aitoff')
    ax2 = fig.add_subplot(212)
    fig.suptitle(title, fontsize=16)

    # Plot visibility cube
    plot_viscube(ax1, args[0])

    # Plot visibility result
    plot_visibility(ax2, args[0])

    # Show figure
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

    # Show visibility cube
    show_viscube()
