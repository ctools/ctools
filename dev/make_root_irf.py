#! /usr/bin/env python
# ==========================================================================
# Generates ROOT IRF test file
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
from ROOT import TFile

# Open input file
filename = '/project-data/cta/performance/prod3b/CTA-Performance-North-E-20170502/CTA-Performance-North-20deg-average-50h_20170502.root'
infile   = TFile(filename)

# Open output file
outfile = TFile('irf.root', 'RECREATE')

# Open histograms
aeff1  = infile.Get('EffectiveAreaEtrueNoTheta2cut_offaxis')
aeff2  = infile.Get('EffectiveAreaEtrueNoTheta2cut')
psf1   = infile.Get('AngRes_offaxis')
psf2   = infile.Get('AngRes80_offaxis')
psf3   = infile.Get('AngularPSF2DEtrue_offaxis')
edisp1 = infile.Get('EestOverEtrue_offaxis')
edisp2 = infile.Get('MigMatrixNoTheta2cut_offaxis')
bkg1   = infile.Get('BGRatePerSqDeg_offaxis')
bkg2   = infile.Get('BGRatePerSqDeg')

# Write histograms
infile.GetList().Write()

# Close output file
outfile.Close()
