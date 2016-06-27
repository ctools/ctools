#! /usr/bin/env python
# ==========================================================================
# Generates Prod3 IRFs from ROOT performance files.
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
import gammalib
import cscripts


# ======================= #
# Set Prod3 DESY database #
# ======================= #
def set_prod3_desy():
    """
    Set Prod3 DESY database

    This function sets the calibration database. By adding dictionaries to
    the db list a number of calibrations can be defined for which IRFs
    should be generated. All calibrations will be added to the calibration
    database.
    """
    # Set database attributes
    #path       = '/project-data/cta/performance/prod3'
    #file       = 'DESY.d20160304.V5.ID0NIM2.prod3-paranalp05-NN.S.3HB1-ND-3.180000s.root'
    path       = 'data'
    file       = 'irf.root'
    psftype    = 'Gauss'
    oversample = 3
    #norm1d     = True   # Will normalise the on-axis 2D IRFs on the 1D IRFs
    norm1d     = False
    bgdinfill  = False
    bgdethres  = 10.0

    # Set database content
    db = [{'inst': 'prod3', 'id': 'South_50h', 'psftype': psftype,
           'oversample': oversample, 'norm1d': norm1d,
           'bgdinfill': bgdinfill, 'bgdethres': bgdethres,
           'path': path, 'file': file}
         ]

    # Return database
    return db


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Get database
    entries = set_prod3_desy()

    # Loop over entries
    for entry in entries:

        # Set optional attributes
        if entry.has_key('rebin'):
            rebin = entry['rebin']
        else:
            rebin = False
        if entry.has_key('eascale'):
            eascale = entry['eascale']
        else:
            eascale = 1.0
        if entry.has_key('bgdscale'):
            bgdscale = entry['bgdscale']
        else:
            bgdscale = 1.0
        if entry.has_key('bgdethres'):
            bgdethres = entry['bgdethres']
        else:
            bgdethres = 1000.0
        if entry.has_key('bgdinfill'):
            bgdinfill = entry['bgdinfill']
        else:
            bgdinfill = False
        if entry.has_key('oversample'):
            oversample = entry['oversample']
        else:
            oversample = 1
        if entry.has_key('norm1d'):
            norm1d = entry['norm1d']
        else:
            norm1d = False

        # Set filename
        filename = entry['path']+'/'+entry['file']

        # Set-up csroot2caldb
        caldb = cscripts.csroot2caldb()
        caldb['infile']        = filename
        caldb['outdir']        = 'prod3'
        caldb['inst']          = entry['inst']
        caldb['id']            = entry['id']
        caldb['version']       = entry['inst']
        caldb['analysis']      = 'DESY'
        caldb['psftype']       = entry['psftype']
        caldb['norm1d']        = norm1d
        caldb['rebin']         = rebin
        caldb['eascale']       = eascale
        caldb['bgdscale']      = bgdscale
        caldb['bgdoversample'] = oversample
        caldb['bgdethres']     = bgdethres
        caldb['bgdinfill']     = bgdinfill
        caldb['logfile']       = 'generate_prod3_irfs.log'

        # Add CALDB entry
        caldb.execute()
