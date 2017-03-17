#! /usr/bin/env python
# ==========================================================================
# Display pointings in observation definition XML file
#
# Copyright (C) 2017 Juergen Knoedlseder
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


# =============================== #
# Extract pointings from XML file #
# =============================== #
def get_pointings(filename):
    """
    Extract pointings from XML file

    Parameters
    ----------
    filename : str
        File name of observation definition XML file

    Returns
    -------
    pnt : list of dict
        Pointings
    """
    # Initialise pointings
    pnt = []

    # Open XML file
    xml = gammalib.GXml(filename)

    # Get observation list
    obs = xml.element('observation_list')

    # Get number of observations
    nobs = obs.elements('observation')

    # Loop over observations
    for i in range(nobs):

        # Get observation
        run = obs.element('observation', i)

        # Get pointing parameter
        npars = run.elements('parameter')
        ra    = None
        dec   = None
        south = False
        for k in range(npars):
            par = run.element('parameter', k)
            if par.attribute('name') == 'Pointing':
                ra  = float(par.attribute('ra'))
                dec = float(par.attribute('dec'))
            if par.attribute('name') == 'Calibration':
                if 'South' in par.attribute('response'):
                    south = True
                else:
                    south = False
        if ra != None:
            p   = gammalib.GSkyDir()
            p.radec_deg(ra, dec)
            entry = {'l': p.l_deg(), 'b': p.b_deg(), 'ra': ra, 'dec': dec, 'south': south}
            pnt.append(entry)

    # Return pointings
    return pnt


# ================ #
# Plot information #
# ================ #
def plot_pointings(pnt, plotfile):
    """
    Plot information

    Parameters
    ----------
    pnt : list of dict
        Pointings
    plotfile : str
        Plot filename
    """
    # Create figure
    plt.figure()

    # Setup figure
    ax = plt.gca()
    ax.cla()
    ax.set_xlim((180, -180))
    ax.set_ylim((-90, 90))

    # Loop over pointings
    for p in pnt:

        # Get longitude
        l = p['l']
        if l > 180.0:
            l = l - 360.0

        # Set color
        if p['south']:
            color='r'
        else:
            color='b'

        # Set circle
        circle = plt.Circle((l, p['b']), 3.0, color=color, fill=False)

        # Add circle
        ax.add_artist(circle)

    # Plot title and labels
    plt.xlabel('Galactic longitude (deg)')
    plt.ylabel('Galactic latitude (deg)')
    plt.title('Pointings')

    # Show plots or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ============== #
# Show pointings #
# ============== #
def show_pointings():
    """
    Show pointings
    """
    # Set usage string
    usage = 'show_pointings.py [-p plotfile] file'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Get pointings
    pnt = get_pointings(args[0])

    # Plot pointings
    plot_pointings(pnt, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show pointings
    show_pointings()
