#!/bin/sh
# =====================================================================
# Run all test and example scripts that come with the package.
#
# Copyright (C) 2011 Jurgen Knodlseder
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
# =====================================================================

#
# Save current working directory
# ==============================
base=$PWD


#
# test
# ====
echo
echo "=====> test"
cd test
./example_models.py
./example_survey.py
cd $base


#
# examples
# ========
echo
echo "=====> examples"
cd examples

# Create local environment
rm -rf gammalib
mkdir -p gammalib/share
mkdir -p gammalib/share/caldb
ln -s $base/caldb gammalib/share/caldb/cta
ln -s $base/models gammalib/share/models
export GAMMALIB=./gammalib

# Run examples
./make_binned_analysis.py
./make_unbinned_analysis.py
#./make_ts_distributions.py
#./make_pull_at_sensitivity_limit.py
cd $base


#
# Signal completion
# =================
echo 
echo "All scripts completed."
