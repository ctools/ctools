#! /usr/bin/env python
# ==========================================================================
# Generates ROOT IRF test file
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
from ROOT import TFile, TH1F, TH2F


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':

    # Set ROOT performance filename
    filename = '/project-data/cta/performance/prod3/Staging.20160518/evndisp/Reference_1_2_180000s_South.root'

    # Open ROOT performance file
    infile = TFile(filename)

    # Get histograms
    ea_true_offaxis     = infile.Get('EffectiveAreaEtrue_offaxis')
    ea_reco_offaxis     = infile.Get('EffectiveArea_offaxis')
    ea_true             = infile.Get('EffectiveAreaEtrue')
    ea_reco             = infile.Get('EffectiveArea')
    angres_offaxis      = infile.Get('AngRes_offaxis')
    angres80_offaxis    = infile.Get('AngRes80_offaxis')
    erecoOetrue_offaxis = infile.Get('EestOverEtrue_offaxis')
    bgd_offaxis         = infile.Get('BGRatePerSqDeg_offaxis')
    bgd                 = infile.Get('BGRatePerSqDeg')
    

    # Open output file
    outfile = TFile('irf.root', 'RECREATE')

    # Write histograms
    ea_true_offaxis.Write()
    ea_reco_offaxis.Write()
    ea_true.Write()
    ea_reco.Write()
    angres_offaxis.Write()
    angres80_offaxis.Write()
    erecoOetrue_offaxis.Write()
    bgd_offaxis.Write()
    bgd.Write()

    # Close file
    outfile.Close()
