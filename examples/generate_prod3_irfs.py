#! /usr/bin/env python
# ==========================================================================
# Generates Prod3 IRFs from ROOT performance files
#
# Copyright (C) 2016 Juergen Knoedlseder
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
import cscripts


# ======================= #
# Set Prod3 DESY database #
# ======================= #
def set_prod3_desy(datadir):
    """
    Set Prod3 DESY database

    This function sets the calibration database. By adding dictionaries to
    the "db" list a number of calibrations can be defined for which IRFs
    should be generated. All calibrations will be added to the calibration
    database.

    Parameters
    ----------
    datadir : str
        Data directory

    Returns
    -------
    db : list of dict
        Database
    """
    # Set database content
    db = [{'inst': 'prod3', 'id': 'South_50h', 'analysis': 'DESY',
           'oversample': 3, 'norm1d': False, 'bgdinfill': False,
           'bgdethres': 10.0, 'path': datadir, 'file': 'irf.root'}]

    # Return database
    return db


# =================== #
# Generate Prod3 IRFs #
# =================== #
def generate_prod3_irfs():
    """
    Generate Prod3 IRFs
    """
    # Get optional argument
    if len(sys.argv) == 2:
        datadir = sys.argv[1]
    else:
        datadir = 'data'

    # Get database
    entries = set_prod3_desy(datadir)

    # Loop over entries
    for entry in entries:

        # Set filename
        filename = entry['path']+'/'+entry['file']

        # Set-up csroot2caldb
        caldb = cscripts.csroot2caldb()
        caldb['infile']        = filename
        caldb['outdir']        = 'prod3'
        caldb['inst']          = entry['inst']
        caldb['id']            = entry['id']
        caldb['version']       = entry['inst']
        caldb['analysis']      = entry.setdefault('analysis', 'DESY')
        caldb['zenith']        = entry.setdefault('zenith', 20.0)
        caldb['azimuth']       = entry.setdefault('azimuth', 0.0)
        caldb['psftype']       = entry.setdefault('psftype', 'Gauss')
        caldb['norm1d']        = entry.setdefault('norm1d', False)
        caldb['rebin']         = entry.setdefault('rebin', False)
        caldb['eascale']       = entry.setdefault('eascale', 1.0)
        caldb['bgdscale']      = entry.setdefault('bgdscale', 1.0)
        caldb['bgdoversample'] = entry.setdefault('oversample', 1)
        caldb['bgdethres']     = entry.setdefault('bgdethres', 1000.0)
        caldb['bgdinfill']     = entry.setdefault('bgdinfill', False)
        caldb['logfile']       = 'generate_prod3_irfs.log'

        # Add CALDB entry
        caldb.execute()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Generate Prod3 IRFs
    generate_prod3_irfs()
