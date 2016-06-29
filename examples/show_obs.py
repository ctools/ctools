#! /usr/bin/env python
# ==========================================================================
# Display summary of observation definition XML file
#
# Required 3rd party modules:
# - matplotlib
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
import ctools
import cscripts
try:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.close()
except:
    print('This script needs the "matplotlib" module')
    sys.exit()


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Print help if wrong number of arguments provided
    if not len(sys.argv) > 1:
        msg = 'Usage: inspect_obs.py <inobs.xml> <ra> <dec>\n'
        msg += 'or\n'
        msg += 'inspect_obs.py <inobs.xml>'
        sys.exit(msg)
        
    # Initialise flag if plots should be shown or saved
    save = False
    if 'save' in sys.argv:
        save = True
        
    # Get input observation XML file
    inobs = gammalib.GFilename(sys.argv[1])
    
    # Assert more than one observation (circumvents having a binned
    # observation without e.g. zenith and azimuth)
    if inobs.is_fits():
        raise RuntimeError('Input observation container should be an XML file')
    
    # Run csobsinfo
    info = cscripts.csobsinfo()
    info['offset'] = False
    
    # Add offset if coordinates are provided
    if len(sys.argv) == 4:
        ra = float(sys.argv[2])
        dec = float(sys.argv[3])
        info['ra'] = ra
        info['dec'] = dec
        info['offset'] = True
    info['inobs'] = inobs.url()
    info['debug'] = True
    info.run()
    
    # Retrieve observation info
    zeniths  = info.zeniths()
    azimuths = info.azimuths()
    offsets  = info.offsets()
    ebounds  = info.ebounds()
    gti      = info.gti()
    
    # Plot zenith angle distribution
    plt.figure()
    zmin = min(zeniths)
    zmax = max(zeniths)
    plt.hist(zeniths, bins=30, range=(zmin, zmax), fc='blue')
    plt.xlabel('Zenith Angle [deg]')
    plt.ylabel('Abundance')
    plt.title('Zenith angle distribution')
    if save:
        plt.savefig('Zenith angle distribution.eps')
    
    # Plot azimuth angle distribution
    plt.figure()
    amin = min(azimuths)
    amax = max(azimuths)
    plt.hist(azimuths, bins=30, range=(amin, amax), fc='blue')
    plt.xlabel('Azimuth Angle (deg)')
    plt.ylabel('Abundance')
    plt.title('Azimuth distribution')
    if save:
        plt.savefig('Azimuth distribution.eps')
    
    # Plot offset if possible
    if info['offset'].boolean():
        plt.figure()
        omin = min(offsets)
        omax = max(offsets)
        plt.hist(offsets, bins=30, range=(omin, omax), fc='blue')
        plt.xlabel('Offset from (RA, DEC)=('+str(ra)+','+str(dec)+') (deg)')
        plt.ylabel('Abundance')
        plt.title('Offset distribution')
        if save:
            plt.savefig('Offset distribution.eps')
    
    # Plot energy thresholds if possible
    if ebounds.size():
        
        emin = []
        emax = []
        for i in range(ebounds.size()):
            emin.append(ebounds.emin(i).log10TeV())
            emax.append(ebounds.emax(i).log10TeV())
        plt.figure()
        plt.hist(emin, bins=80, range=(-1.0, 2.0), fc='red', label='emin')
        plt.hist(emax, bins=80, range=(-1.0, 2.0), fc='blue', label='emax')  
        plt.xlabel('Energy threshold (log10 (E/TeV))')
        plt.ylabel('Abundance')
        plt.legend(loc='upper left')
        plt.title('Energy threshold')
        if save:
            plt.savefig('Energy threshold.eps')
    
    # Plot observation point in time wrt zenith angle
    plt.figure()
    times = []
    for i in range(gti.size()):
        tmean = gti.tstart(i) + 0.5*(gti.tstart(i)-gti.tstart(i))
        times.append(tmean.mjd())
    plt.plot(times, zeniths, 'o', lw=2.0, color='black')
    plt.xlabel('Time (MJD)')
    plt.ylabel('Zenith Angle (deg)')
    plt.title('Observation time')
    if save:
        plt.savefig('Observation time.eps')

    # Display plots if possible
    if not save:
        
        # Show plots
        plt.show()
        