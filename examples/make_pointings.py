#! /usr/bin/env python
# ==========================================================================
# This script generated pointing patterns for CTA observations.
#
# Usage:
#   ./make_pointings.py
#
# ==========================================================================
#
# Copyright (C) 2015 Juergen Knoedlseder
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
import gammalib
import sys
import math


# ============================================================ #
# Get positions of pointings on sky for an interleaved pattern #
# ============================================================ #
def get_positions(xmin, xmax, ymin, ymax, step):
    """
    Get positions for a patch of interleaved pointings on the sky.
    The function does respect the latitude dependence of the
    pointing distance. The step parameter is for a latitude of
    zero.
    """
    # Initialise positions
    positions = []

    # Determine ystep and number of rows
    ystep   = step * math.sqrt(0.75)
    ny      = int((ymax-ymin)/ystep+1.5)
    rescale = (ymax-ymin)/(ystep*(ny-1))
    ystep  *= rescale
    print("Number of rows ....: %d" % (ny))

    # Loop over rows
    y = ymin
    for row in range(ny):

        # Determine xstep and number of pointings in row
        xstep = step * math.cos(gammalib.deg2rad * y)
        nx    = int((xmax-xmin)/xstep+1.5)
        rescale = (xmax-xmin)/(xstep*(nx-1))
        xstep  *= rescale

        # Set x offset
        if row % 2 == 0:
            x = xmin
        else:
            x = xmin + 0.5 * xstep
        #print(row, x)

        # Append pointings
        for pnt in range(nx):
            positions.append({'x': x, 'y': y})
            x += xstep
            if x >= 360.0:
                x -= 360.0

        # Increment y
        y += ystep

    # Return positions
    return positions


# ====================== #
# Setup patch on the sky #
# ====================== #
def set_patch(lmin=-30.0, lmax=30.0, bmin=-1.5, bmax=+1.5, \
              separation=3.0, hours=100.0, \
              site=None, lst=True, autodec=0.0):
    """
    Setup pointing patch on the sky.
    """
    # Initialise observation definition
    obsdef = []

    # Get pointing positions
    pointings = get_positions(lmin, lmax, bmin, bmax, separation)

    # Determine number of pointings and pointing duration
    n_pnt    = len(pointings)
    duration = hours*3600.0/float(n_pnt)

    # Observing time collectors
    exposure_south = 0.0
    exposure_north = 0.0

    # Set observations
    for pnt in pointings:

        # Set positions and duration
        obs = {'lon': pnt['x'], 'lat': pnt['y'], 'duration': duration}

        # Optionally add-in site dependent IRF
        if site != None:

            # Handle automatic site switching
            if site == "Automatic":
                pos = gammalib.GSkyDir()
                pos.lb_deg(pnt['x'], pnt['y'])
                dec = pos.dec_deg()
                if (dec >= autodec):
                    obs_site = "North"
                else:
                    obs_site = "South"
            else:
                obs_site = site

            # Set site IRF
            if obs_site == "North":
                caldb = "tenerife"
                if lst:
                    irf = "DESY20140105_50h"
                else:
                    irf = "DESY20140105_50h_noLST"
                exposure_north += duration
            else:
                caldb = "aar"
                if lst:
                    irf = "DESY20140105_50h"
                else:
                    irf = "DESY20140105_50h_noLST"
                exposure_south += duration

            # Set IRF information
            obs['caldb'] = caldb
            obs['irf']   = irf

        # Append observation
        obsdef.append(obs)

    # Dump statistics
    print("Number of pointings: %d (%.2f sec)" % (n_pnt,duration))
    print("South array .......: %.2f hours (%.2f%%)" % \
          (exposure_south/3600.0, \
           100.0*exposure_south/(exposure_south+exposure_north)))
    print("North array .......: %.2f hours (%.2f%%)" % \
          (exposure_north/3600.0, \
           100.0*exposure_north/(exposure_south+exposure_north)))

    # Return observation definition
    return obsdef


# =========================== #
# Setup Galactic plane survey #
# =========================== #
def set_gps(separation=3.0, bmin=-1.3, bmax=1.3, lst=True):
    """
    Setup Galactic plane survey.
    """
    # Initialise observation definition
    obsdef = []

    # Add inner region South
    obsdef.extend(set_patch(lmin=-60.0, lmax=60.0, bmin=bmin, bmax=bmax, \
                            separation=separation, hours=780, \
                            site="South", lst=lst))

    # Add Vela & Carina region
    obsdef.extend(set_patch(lmin=240.0, lmax=300.0, bmin=bmin, bmax=bmax, \
                            separation=separation, hours=180, \
                            site="South", lst=lst))

    # Add 210-240 region
    obsdef.extend(set_patch(lmin=210.0, lmax=240.0, bmin=bmin, bmax=bmax, \
                            separation=separation, hours=60, \
                            site="South", lst=lst))

    # Add Cygnus, Perseus
    obsdef.extend(set_patch(lmin=60.0, lmax=150.0, bmin=bmin, bmax=bmax, \
                            separation=separation, hours=450, \
                            site="North", lst=lst))

    # Add Anticentre
    obsdef.extend(set_patch(lmin=150.0, lmax=210.0, bmin=bmin, bmax=bmax, \
                            separation=separation, hours=150, \
                            site="North", lst=lst))

    # Return observation definition
    return obsdef


# ========================== #
# Setup Extragalactic survey #
# ========================== #
def set_extgal(separation=3.0, lst=True):
    """
    Setup Extragalactic survey.
    """
    # Initialise observation definition
    obsdef = []

    # Set patch
    obsdef.extend(set_patch(lmin=-90.0, lmax=90.0, bmin=+5.0, bmax=+90.0, \
                            separation=separation, hours=1000, \
                            site="Automatic", lst=lst, autodec=25.0))

    # Return observation definition
    return obsdef


# ========================= #
# Setup Galactic centre KSP #
# ========================= #
def set_gc(lst=True):
    """
    Setup Galactic centre KSP.
    """
    # Initialise observation definition
    obsdef = []

    # Central wobble
    obsdef.extend(set_patch(lmin=-0.5, lmax=0.5, bmin=-0.5, bmax=0.5, \
                            separation=0.1, hours=525, \
                            site="South", lst=lst))

    # Extended region
    obsdef.extend(set_patch(lmin=-10.0, lmax=10.0, bmin=-10.0, bmax=10.0, \
                            separation=1.5, hours=300, \
                            site="South", lst=lst))

    # Return observation definition
    return obsdef


# ======================================= #
# Write observation definition dictionary #
# ======================================= #
def write_obsdef(filename, obsdef):
    """
    Write observation definition file.
    """
    # Open file
    file = open(filename, 'w')

    # Write header
    file.write("ra,dec,duration,caldb,irf\n")

    # Loop over pointings
    for obs in obsdef:

        # If we have lon,lat then convert into RA,Dec
        if "lon" in obs and "lat" in obs:
            lon = obs["lon"]
            lat = obs["lat"]
            dir = gammalib.GSkyDir()
            dir.lb_deg(lon,lat)
            ra  = dir.ra_deg()
            dec = dir.dec_deg()
        else:
            ra  = obs["ra"]
            dec = obs["dec"]

        # Write information
        file.write("%8.4f,%8.4f,%.4f,%s,%s\n" % \
                   (ra, dec, obs['duration'], obs['caldb'], obs['irf']))

    # Close file
    file.close()

    # Return
    return
    

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    CTA pointing pattern generation.
    """
    # Initialise flags
    need_help = False
    
    # Test for command line arguments
    #print(sys.argv[0])
    if (len(sys.argv) > 1):
        if sys.argv[1] == "-h":
            need_help = True
        else:
            obsname = sys.argv[1]
    else:
        need_help = True

    # Print help if needed and exit
    if need_help:
        print("Usage: make_pointing.py [OPTIONS]")
        print("     -h       Display this usage message")
        print("     gps      Galactic plane scan (2 row scheme)")
        print("     gps3     Galactic plane scan (3 row scheme)")
        print("     extgal   Extragalactic survey")
        print("     gc       Galactic centre")
        sys.exit()

    # Galactic plane scan
    if obsname == "gps":
        obsdef = set_gps(lst=True)
        write_obsdef("gps.dat", obsdef)
    elif obsname == "gps3":
        obsdef = set_gps(separation=1.5, lst=True)
        write_obsdef("gps3.dat", obsdef)

    # Extragalactic survey
    elif obsname == "extgal":
        obsdef = set_extgal(lst=True)
        write_obsdef("extgal.dat", obsdef)

    # Galactic centre
    elif obsname == "gc":
        obsdef = set_gc(lst=True)
        write_obsdef("gc.dat", obsdef)

    # Invalid pattern
    else:
        print("Unknown option \""+obsname+"\"")
        
