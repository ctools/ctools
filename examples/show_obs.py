#! /usr/bin/env python
# ==========================================================================
# Display summary of observation definition XML file
#
# Copyright (C) 2016 Michael Mayer
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


# ==================== #
# Run csobsinfo script #
# ==================== #
def run_csobsinfo(filename, ra, dec, debug=True):
    """
    Run csobsinfo script

    Parameters
    ----------
    filename : str
        File name of observation definition XML file
    ra : str
        Target Right Ascension (deg)
    dec : str
        Target declination of pointing (deg)
    debug : bool, optional
        Switch on debugging in csobsinfo run

    Returns
    -------
    info : `~cscripts.csobsinfo`
        csobsinfo instance
    """
    # Get observation definition XML filename
    obsdef = gammalib.GFilename(filename)

    # Setup csobsinfo script
    info = cscripts.csobsinfo()
    info['inobs'] = obsdef.url()
    info['debug'] = debug

    # Set offset if the "ra" and "dec" arguments are not empty strings
    if ra != '' and dec != '':
        info['ra']     = float(ra)
        info['dec']    = float(dec)
        info['offset'] = True
    else:
        info['offset'] = False

    # Run csobsinfo
    info.run()

    # Return
    return info


# ================ #
# Plot information #
# ================ #
def plot_information(info, ra, dec, plotfile):
    """
    Plot information

    Parameters
    ----------
    info : `~cscripts.csobsinfo`
        csobsinfo instance
    ra : float
        Right Ascension (deg)
    dec : float
        Declination (deg)
    plotfile : str
        Plot filename
    """
    # Retrieve observation information
    zeniths  = info.zeniths()
    azimuths = info.azimuths()
    offsets  = info.offsets()
    ebounds  = info.ebounds()
    gti      = info.gti()

    # Create figure with subplots. In case that offset and energy boundary
    # information is available a column is added to the subplots
    nrows = 2
    ncols = 2
    if info['offset'].boolean() and ebounds.size() > 0:
        ncols += 1
    iplot = 1
    plt.figure()
    plt.subplot(nrows, ncols, iplot)

    # Plot zenith angle histogram
    zmin = min(zeniths)
    zmax = max(zeniths)
    plt.hist(zeniths, bins=30, range=(zmin, zmax), fc='blue')
    plt.xlabel('Zenith Angle (deg)')
    plt.ylabel('Frequency')
    plt.title('Zenith angle distribution')

    # Plot azimuth angle histogram
    iplot += 1
    plt.subplot(nrows, ncols, iplot)
    amin = min(azimuths)
    amax = max(azimuths)
    plt.hist(azimuths, bins=30, range=(amin, amax), fc='blue')
    plt.xlabel('Azimuth Angle (deg)')
    plt.ylabel('Frequency')
    plt.title('Azimuth distribution')

    # If available, plot offset from target histogram
    if info['offset'].boolean():
        iplot += 1
        plt.subplot(nrows, ncols, iplot)
        omin = min(offsets)
        omax = max(offsets)
        plt.hist(offsets, bins=30, range=(omin, omax), fc='blue')
        plt.xlabel('Offset from (RA, DEC)=('+str(ra)+','+str(dec)+') (deg)')
        plt.ylabel('Frequency')
        plt.title('Offset distribution')

    # If available, plot energy thresholds histogram
    if ebounds.size() > 0:
        emin = []
        emax = []
        for i in range(ebounds.size()):
            emin.append(ebounds.emin(i).log10TeV())
            emax.append(ebounds.emax(i).log10TeV())
        iplot += 1
        plt.subplot(nrows, ncols, iplot)
        plt.hist(emin, bins=80, range=(-1.0, 2.0), fc='red', label='emin')
        plt.hist(emax, bins=80, range=(-1.0, 2.0), fc='blue', label='emax')  
        plt.xlabel('Energy threshold (log10 (E/TeV))')
        plt.ylabel('Frequency')
        plt.legend(loc='upper left')
        plt.title('Energy threshold')

    # Plot zenith angle versus mean observing time for all observations
    iplot += 1
    plt.subplot(nrows, ncols, iplot)
    times = []
    for i in range(gti.size()):
        tmean = gti.tstart(i) + 0.5*(gti.tstart(i)-gti.tstart(i))
        times.append(tmean.mjd())
    plt.plot(times, zeniths, 'o', lw=2.0, color='black')
    plt.xlabel('Time (MJD)')
    plt.ylabel('Zenith Angle (deg)')
    plt.title('Observation time')

    # Show plots or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ================ #
# Show observation #
# ================ #
def show_obs():
    """
    Show observation
    """
    # Set usage string
    usage = 'show_obs.py [-ra ra] [-dec dec] [-p plotfile] file'

    # Set default options
    options = [{'option': '-ra',  'value': ''},
               {'option': '-dec', 'value': ''},
               {'option': '-p',   'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    ra       = options[0]['value']
    dec      = options[1]['value']
    plotfile = options[2]['value']

    # Run csobsinfo
    info = run_csobsinfo(args[0], ra, dec)

    # Plot information
    plot_information(info, ra, dec, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show observation
    show_obs()
