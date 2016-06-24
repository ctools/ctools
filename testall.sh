#!/bin/sh
# =====================================================================
# Run some test scripts that come with the package
#
# Copyright (C) 2011-2016 Juergen Knoedlseder
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
# test/dev
# ========
echo
echo "=====> test/dev"
cd test/dev

# Create local environment
rm -rf ctools
mkdir -p ctools/share/caldb
ln -s $base/caldb ctools/share/caldb/cta
ln -s $base/models ctools/share/models
export CTOOLS=$PWD/ctools

# Run checkers
./check_models.py

# Run makers
./make_survey.py

# Remove local environment
rm -rf ctools

# Step back to base directory
cd $base


#
# Signal completion
# =================
echo 
echo "All scripts completed."
