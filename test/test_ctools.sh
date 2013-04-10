#! /bin/sh
# ==========================================================================
# This script tests all ctools
#
# Copyright (C) 2011-2012 Juergen Knoedlseder
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
# Get echo command (to circumvent built-in echo on some systems)
#
ECHO=`which echo`

#
# Print Header
#
$ECHO "***************"
$ECHO "* Test ctools *"
$ECHO "***************"

#
# Define executables
# ==================
ctbin=../src/ctbin/ctbin
ctlike=../src/ctlike/ctlike
ctmodel=../src/ctmodel/ctmodel
ctobssim=../src/ctobssim/ctobssim
ctselect=../src/ctselect/ctselect
ctskymap=../src/ctskymap/ctskymap


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
$ECHO -n "Test ctobssim: "
$ctobssim infile="data/crab.xml" \
          outfile="events.fits" \
          caldb="irf" \
          irf="cta_dummy_irf" \
          ra=83.63 \
          dec=22.01 \
          rad=10.0 \
          tmin=0.0 \
          tmax=1800.0 \
          emin=0.1 \
          emax=100.0
$ECHO -n "."
if [ -s "events.fits" ]
then
  $ECHO -n "."
else
  echo " events.fits file is not valid"
  exit 1
fi
$ECHO " ok"


#
# Test ctskymap
# =============
$ECHO -n "Test ctskymap: "
$ctskymap evfile="events.fits" \
          outfile="skymap.fits" \
          emin=0.1 \
          emax=100.0 \
          nxpix=200 \
          nypix=200 \
          binsz=0.02 \
          coordsys="CEL" \
          xref=83.63 \
          yref=22.01 \
          proj="CAR"
$ECHO -n "."
if [ -s "skymap.fits" ]
then
  $ECHO -n "."
else
  $ECHO " skymap.fits file is not valid"
  exit 1
fi
$ECHO " ok"


#
# Test ctbin
# ==========
$ECHO -n "Test ctbin: "
$ctbin evfile="events.fits" \
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
$ECHO -n "."
if [ -s "cntmap.fits" ]
then
  $ECHO -n "."
else
  $ECHO " cntmap.fits file is not valid"
  exit 1
fi
$ECHO " ok"


#
# Test ctmodel
# ============
$ECHO -n "Test ctmodel: "
#
# Run 1
$ctmodel infile="cntmap.fits" \
         outfile="modmap1.fits" \
         caldb="irf" \
         irf="cta_dummy_irf" \
         srcmdl="data/crab.xml"
$ECHO -n "."
if [ -s "modmap1.fits" ]
then
  $ECHO -n "."
else
  $ECHO " modmap1.fits file is not valid"
  exit 1
fi
#
# Run 2
$ctmodel infile="NONE" \
         outfile="modmap2.fits" \
         caldb="irf" \
         irf="cta_dummy_irf" \
         srcmdl="data/crab.xml" \
         ra=83.63 \
         dec=22.01 \
         tmin=0.0 \
         tmax=1800.0 \
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
$ECHO -n "."
if [ -s "modmap2.fits" ]
then
  $ECHO -n "."
else
  $ECHO " modmap2.fits file is not valid"
  exit 1
fi
$ECHO " ok"


#
# Test ctselect
# =============
$ECHO -n "Test ctselect: "
$ctselect infile="events.fits" \
          outfile="selected_events.fits" \
          ra=83.63 \
          dec=22.01 \
          rad=10.0 \
          tmin=0.0 \
          tmax=1800.0 \
          emin=0.1 \
          emax=100.0
$ECHO -n "."
if [ -s "selected_events.fits" ]
then
  $ECHO -n "."
else
  $ECHO " selected_events.fits file is not valid"
  exit 1
fi
$ECHO " ok"


#
# Test ctlike
# ===========
$ECHO -n "Test ctlike: "
$ctlike infile="cntmap.fits" \
        srcmdl="data/crab.xml" \
        outmdl="results_binned.xml" \
        caldb="irf" \
        irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "results_binned.xml" ]
then
  $ECHO -n "."
else
  $ECHO " results_binned.xml file is not valid"
  exit 1
fi
$ctlike infile="selected_events.fits" \
        srcmdl="data/crab.xml" \
        outmdl="results_unbinned.xml" \
        caldb="irf" \
        irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "results_unbinned.xml" ]
then
  $ECHO -n "."
else
  $ECHO " results_unbinned.xml file is not valid"
  exit 1
fi
$ECHO " ok"
