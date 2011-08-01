#! /bin/bash
# ==========================================================================
# This script tests all cscripts that are shipped with ctatools
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
#
# ==========================================================================

#
# Print Header
#
echo ""
echo "*****************"
echo "* Test cscripts *"
echo "*****************"


#
# Remove any existing result files
# ================================
rm -rf *.fits *.log *.xml *.dat


#
# Creates pfiles directory
# ========================
mkdir -p pfiles
export PFILES=pfiles


#
# Test cssens
# ===========
echo -n "Test cssens: "
rm -rf cssens
ln -s ../scripts/cssens.py cssens
chmod +x ../scripts/cssens.py
cssens duration=3600.0 \
       caldb="irf" \
       irf="kb_E_50h_v3" \
       type="point" \
       offset="0.0" \
       bkg="data/bkg_kb_E_50h_v3.txt" \
       rad=1.0
echo -n "."
if [ -s "sensitivity.dat" ]
then
  echo -n "."
else
  echo " sensitivity.dat file is not valid"
  exit 1
fi
echo " ok"
