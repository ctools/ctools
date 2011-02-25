#!/bin/sh
#
# Run all test and example scripts that come with the package
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
base=$PWD


# test
# ====
echo
echo "=====> test"
cd test
./test_coverage.py "data/crab.xml" 10 pull.dat
./example_binned.py
./exmaple_models.py
./example_survey.py
./example_unbinned.py
cd $base


# Signal completion
# =================
echo 
echo "All scripts completed."
