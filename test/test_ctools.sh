#! /bin/bash
# ==========================================================================
# This script tests all ctools
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
echo "***************"
echo "* Test ctools *"
echo "***************"


#
# Remove any existing result files
# ================================
rm -rf *.fits *.log *.xml


#
# Creates pfiles directory
# ========================
rm -rf pfiles
mkdir -p pfiles
cp -r ../src/*/*.par pfiles/
export PFILES=pfiles


#
# Test ctobssim
# =============
echo -n "Test ctobssim: "
ctobssim infile="data/crab.xml" \
         outfile="events.fits" \
         caldb="irf" \
         irf="kb_E_50h_v3" \
         ra=83.63 \
         dec=22.01 \
         rad=10.0 \
         tmin=0.0 \
         tmax=1800.0 \
         emin=0.1 \
         emax=100.0
echo -n "."
if [ -s "events.fits" ]
then
  echo -n "."
else
  echo " events.fits file is not valid"
  exit 1
fi
echo " ok"


#
# Test ctbin
# ==========
echo -n "Test ctbin: "
ctbin evfile="events.fits" \
      outfile="cntmap.fits" \
      emin=0.1 \
      emax=100.0 \
      enumbins=20 \
      nxpix=200 \
      nypix=200 \
      binsz=0.02 \
      coordsys="CEL" \
      xref=83.63 \
      yref=22.01 \
      proj="CAR"
echo -n "."
if [ -s "cntmap.fits" ]
then
  echo -n "."
else
  echo " cntmap.fits file is not valid"
  exit 1
fi
echo " ok"


#
# Test ctselect
# =============
echo -n "Test ctselect: "
ctselect infile="events.fits" \
         outfile="selected_events.fits" \
         ra=83.63 \
         dec=22.01 \
         rad=10.0 \
         tmin=0.0 \
         tmax=1800.0 \
         emin=0.1 \
         emax=100.0
echo -n "."
if [ -s "selected_events.fits" ]
then
  echo -n "."
else
  echo " selected_events.fits file is not valid"
  exit 1
fi
echo " ok"


#
# Test ctlike
# ===========
echo -n "Test ctlike: "
ctlike cntmap="cntmap.fits" \
       srcmdl="data/crab.xml" \
       outmdl="results_binned.xml" \
       method="BINNED" \
       caldb="irf" \
       irf="kb_E_50h_v3"
echo -n "."
if [ -s "results_binned.xml" ]
then
  echo -n "."
else
  echo " results_binned.xml file is not valid"
  exit 1
fi
ctlike evfile="selected_events.fits" \
       srcmdl="data/crab.xml" \
       outmdl="results_unbinned.xml" \
       method="UNBINNED" \
       caldb="irf" \
       irf="kb_E_50h_v3"
echo -n "."
if [ -s "results_unbinned.xml" ]
then
  echo -n "."
else
  echo " results_unbinned.xml file is not valid"
  exit 1
fi
echo " ok"
