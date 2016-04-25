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
    Set Prod3 DESY database.
    """
    # Set database attributes
    path    = "/project-data/cta/performance/prod3"
    psftype = "Gauss"

    # Set database content
    db = [{'inst': "prod3", 'id': "South_50h", 'psftype': psftype, 'path': path,
           'file': "DESY.d20160304.V5.ID0NIM2.prod3-paranalp05-NN.S.3HB1-ND-3.180000s.root"}
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

        # Set filename
        filename = entry['path']+"/"+entry['file']

        # Set-up csroot2caldb
        caldb = cscripts.csroot2caldb()
        caldb["infile"]   = filename
        caldb["outdir"]   = "prod3"
        caldb["inst"]     = entry['inst']
        caldb["id"]       = entry['id']
        caldb["version"]  = entry['inst']
        caldb["analysis"] = "DESY"
        caldb["psftype"]  = entry['psftype']
        caldb["rebin"]    = rebin
        caldb["eascale"]  = eascale
        caldb["bgdscale"] = bgdscale
        caldb["logfile"]  = "generate_prod3_irfs.log"

        # Add CALDB entry
        caldb.execute()
