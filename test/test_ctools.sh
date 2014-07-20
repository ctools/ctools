#! /bin/sh
# ==========================================================================
# This script tests all ctools
#
# Copyright (C) 2011-2014 Juergen Knoedlseder
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
ctexpcube=../src/ctexpcube/ctexpcube
ctcubemask=../src/ctcubemask/ctcubemask


#
# Remove any existing result files
# ================================
rm -rf *.fits *.xml


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
#
# Run 1
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
#
# Run 2
$ctobssim infile="data/model_background_cube.xml" \
          outfile="events2.fits" \
          caldb="irf" \
          irf="cta_dummy_irf" \
          ra=84.17263 \
          dec=22.01444 \
          rad=1.0 \
          tmin=0.0 \
          tmax=1.0 \
          emin=0.3 \
          emax=8.0 chatter=4
$ECHO -n "."
if [ -s "events2.fits" ]
then
  $ECHO -n "."
else
  echo " events2.fits file is not valid"
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
#
# Run 1
$ctbin evfile="events.fits" \
       outfile="cntmap.fits" \
       emin=0.1 \
       emax=100.0 \
       enumbins=20 \
       ebinalg="LOG" \
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
#
# Run 2
$ctbin evfile="events2.fits" \
       outfile="cntmap2.fits" \
       emin=0.3 \
       emax=8.0 \
       enumbins=10 \
       ebinalg="LOG" \
       nxpix=50 \
       nypix=50 \
       binsz=0.02 \
       coordsys="CEL" \
       xref=84.17263 \
       yref=22.01444 \
       proj="CAR"
$ECHO -n "."
if [ -s "cntmap2.fits" ]
then
  $ECHO -n "."
else
  $ECHO " cntmap2.fits file is not valid"
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
# Run 3
$ctmodel infile="NONE" \
         outfile="modmap3.fits" \
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
if [ -s "modmap3.fits" ]
then
  $ECHO -n "."
else
  $ECHO " modmap3.fits file is not valid"
  exit 1
fi
#
# Run 2
$ctmodel infile="cntmap2.fits" \
         outfile="modmap2.fits" \
         caldb="irf" \
         irf="cta_dummy_irf" \
         srcmdl="data/model_background_cube.xml"
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
#
# Run 1
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
#
# Run 2
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
#
# Run 3
$ctlike infile="events2.fits" \
        srcmdl="data/model_background_cube.xml" \
        outmdl="results_unbinned_background_cube.xml" \
        caldb="irf" \
        irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "results_unbinned_background_cube.xml" ]
then
  $ECHO -n "."
else
  $ECHO " results_unbinned_background_cube.xml file is not valid"
  exit 1
fi
#
# Run 4
$ctlike infile="cntmap2.fits" \
        srcmdl="data/model_background_cube.xml" \
        outmdl="results_binned_background_cube.xml" \
        caldb="irf" \
        irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "results_binned_background_cube.xml" ]
then
  $ECHO -n "."
else
  $ECHO " results_binned_background_cube.xml file is not valid"
  exit 1
fi
$ECHO " ok"


#
# Test ctexpcube
# ==============
$ECHO -n "Test ctexpcube: "
#
# Run 1
$ctexpcube infile="data/crab_events.fits" \
           cntmap="NONE" \
           outfile="expcube1.fits" \
           caldb="irf" \
           irf="cta_dummy_irf" \
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
if [ -s "expcube1.fits" ]
then
  $ECHO -n "."
else
  $ECHO " expcube1.fits file is not valid"
  exit 1
fi
#
# Run 2
$ctexpcube infile="data/crab_events.fits" \
           cntmap="data/crab_cntmap.fits" \
           outfile="expcube2.fits" \
           caldb="irf" \
           irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "expcube2.fits" ]
then
  $ECHO -n "."
else
  $ECHO " expcube2.fits file is not valid"
  exit 1
fi
$ECHO " ok"


#
# Test ctcubemask
# ===============
$ECHO -n "Test ctcubemask: "
$ctcubemask infile="data/crab_events.fits" \
            regfile="data/exclusion.reg" \
            outfile="filtered_cube.fits" \
            ra=83.63 \
            dec=22.01 \
            rad=2.0 \
            emin=0.1 \
            emax=100.0
$ECHO -n "."
if [ -s "filtered_cube.fits" ]
then
  $ECHO -n "."
else
  $ECHO " filtered_cube.fits file is not valid"
  exit 1
fi
$ECHO " ok"
