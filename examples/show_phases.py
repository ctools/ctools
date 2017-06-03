#! /usr/bin/env python
# ==========================================================================
# Display phase information from event list
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


# =========== #
# Plot phases #
# =========== #
def plot_phases(filename, plotfile, nphases=40):
    """
    Plot phases

    Parameters
    ----------
    filename : str
        Name of lightcurve FITS file
    plotfile : str
        Plot file name
    """
    # Read observation container. If an exception occurs then try loading the
    # file as an event list and build an observation container.
    try:
        obs = gammalib.GObservations(filename)
    except:
        obs = gammalib.GObservations()
        obs.append(gammalib.GCTAObservation(filename))

    # Initialise phases
    phases = []

    # Loop over all observations
    for run in obs:

        # Get events
        events = run.events()

        # Loop over all events
        for event in events:

            # Get phase
            phase = event.phase()

            # Collect phase
            phases.append(phase)
            phases.append(phase+1.0)

    # Plot phase
    plt.figure()
    plt.hist(phases,bins=nphases,range=(0.0,2.0),facecolor='r')
    plt.xlabel('Phase')
    plt.ylabel('Events')

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# =========== #
# Show phases #
# =========== #
def show_phases():
    """
    Show phases
    """
    # Set usage string
    usage = 'show_phases.py [-p plotfile] [file]'

    # Set default options
    options = [{'option': '-n', 'value': '20'},
               {'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    nphases  = int(options[0]['value'])*2
    plotfile = options[1]['value']

    # Plot phases
    plot_phases(args[0], plotfile, nphases=nphases)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show phases
    show_phases()
