#! /usr/bin/env python
# ==========================================================================
# This script simulates pointing patterns for CTA surveys.
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


# ======================= #
# Setup double row scheme #
# ======================= #
def set_double_row(lmin=-30.0, lmax=30.0, separation=3.0, hours=100.0, site=None):
    """
    Setup double row scheme.
    """
    # Initialise observation definition
    obsdef = []

    # Compute number of pointings in row
    n_pnt_per_row = int((lmax-lmin)/separation+0.5)
    n_pnt         = 2*n_pnt_per_row

    # Compute offset from Galactic plane (degrees)
    offset = separation * 0.5 * math.sqrt(0.75)

    # Compute duration per pointing (seconds)
    duration = hours*3600.0/float(n_pnt)

    # Set upper row
    lon = lmin
    for i in range(n_pnt_per_row):
        obs = {'lon': lon, 'lat': offset, 'duration': duration}
        obsdef.append(obs)
        lon += separation

    # Set lower row
    lon = lmin + 0.5*separation
    for i in range(n_pnt_per_row):
        obs = {'lon': lon, 'lat': -offset, 'duration': duration}
        obsdef.append(obs)
        lon += separation

    # Optionally add-in site dependent IRF
    if site != None:
        if site == "North":
            caldb = "tenerife"
            irf   = "DESY20140105_50h"
        else:
            caldb = "aar"
            irf   = "DESY20140105_50h"
        for obs in obsdef:
            obs['caldb'] = caldb
            obs['irf']   = irf

    # Return observation definition
    return obsdef


# =========================== #
# Setup Galactic plane survey #
# =========================== #
def set_gps():
    """
    Setup Galactic plane survey.
    """
    # Initialise observation definition
    obsdef = []

    # Add inner region South
    obsdef.extend(set_double_row(lmin=-60.0, lmax=60.0, hours=780, site="South"))

    # Add Vela & Carina region
    obsdef.extend(set_double_row(lmin=240.0, lmax=300.0, hours=180, site="South"))

    # Add 210-240 region
    obsdef.extend(set_double_row(lmin=210.0, lmax=240.0, hours=60, site="South"))

    # Add Cygnus, Perseus
    obsdef.extend(set_double_row(lmin=60.0, lmax=150.0, hours=450, site="North"))

    # Add Anticentre
    obsdef.extend(set_double_row(lmin=150.0, lmax=210.0, hours=150, site="North"))

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
    CTA survey simulation.
    """
    # Initialise flags
    need_help = False
    
    # Test for command line arguments
    print(sys.argv[0])
    if (len(sys.argv) > 1):
        if sys.argv[1] == "-h":
            need_help = True
        else:
            need_help = True

    # Print help if needed and exit
    if need_help:
        print("Usage: make_pointing.py [OPTIONS]")
        print("     -h       Display this usage message")
        sys.exit()

    # Test
    obsdef = set_gps()

    # Write observation definition file
    write_obsdef("gps.dat", obsdef)
