#! /usr/bin/env python
# ==========================================================================
# Generates test data
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
import gammalib
import ctools


# ================= #
# Create event list #
# ================= #
def create_event_list():
    """
    Create event list
    """
    # Simulate events
    sim = ctools.ctobssim()
    sim['inmodel']   = 'test/data/crab.xml'
    sim['outevents'] = 'test/data/crab_events.fits'
    sim['caldb']     = 'prod2'
    sim['irf']       = 'South_0.5h'
    sim['edisp']     = False
    sim['ra']        = 83.63
    sim['dec']       = 22.51
    sim['rad']       = 2.0
    sim['tmin']      = '2020-01-01T00:00:00'
    sim['tmax']      = '2020-01-01T00:05:00'
    sim['emin']      = 1.0
    sim['emax']      = 100.0
    sim.execute()

    # Return
    return


# ================== #
# Create counts cube #
# ================== #
def create_count_cube():
    """
    Create counts cube
    """
    # Allocate ctbin application and set parameters
    bin = ctools.ctbin()
    bin['inobs']   = 'test/data/crab_events.fits'
    bin['outcube'] = 'test/data/crab_cntmap.fits'
    bin['ebinalg']  = 'LOG'
    bin['emin']     = 1.0
    bin['emax']     = 100.0
    bin['enumbins'] = 5
    bin['usepnt']   = True # Use pointing for map centre
    bin['nxpix']    = 20
    bin['nypix']    = 20
    bin['binsz']    = 0.2
    bin['coordsys'] = 'CEL'
    bin['proj']     = 'TAN'
    bin.execute()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run tests
    create_event_list()
    create_count_cube()

