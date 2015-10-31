#! /bin/sh
# ==========================================================================
# This script tests all ctools
#
# Copyright (C) 2011-2015 Juergen Knoedlseder
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
cttsmap=../src/cttsmap/cttsmap
ctmodel=../src/ctmodel/ctmodel
ctobssim=../src/ctobssim/ctobssim
ctselect=../src/ctselect/ctselect
ctskymap=../src/ctskymap/ctskymap
ctexpcube=../src/ctexpcube/ctexpcube
ctpsfcube=../src/ctpsfcube/ctpsfcube
ctbkgcube=../src/ctbkgcube/ctbkgcube
ctcubemask=../src/ctcubemask/ctcubemask
ctbutterfly=../src/ctbutterfly/ctbutterfly
ctulimit=../src/ctulimit/ctulimit
cterror=../src/cterror/cterror


#
# Remove any existing result files
# ================================
rm -rf *.fits *.xml
rm -rf ulimit.dat butterfly.txt


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
$ctobssim inmodel="data/crab.xml" \
          outevents="events.fits" \
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
  echo " events.fits file is not found"
  exit 1
fi
$ECHO " ok"


#
# Test ctskymap
# =============
$ECHO -n "Test ctskymap: "
$ctskymap inobs="events.fits" \
          outmap="skymap.fits" \
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
  $ECHO " skymap.fits file is not found"
  exit 1
fi
$ECHO " ok"


#
# Test ctbin
# ==========
$ECHO -n "Test ctbin: "
#
# Run 1
$ctbin inobs="events.fits" \
       outcube="cntmap1.fits" \
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
if [ -s "cntmap1.fits" ]
then
  $ECHO -n "."
else
  $ECHO " cntmap1.fits file is not found"
  exit 1
fi
#
# Run 2
$ctbin inobs="events.fits" \
       outcube="cntmap2.fits" \
       emin=0.3 \
       emax=8.0 \
       enumbins=7 \
       ebinalg="LOG" \
       nxpix=30 \
       nypix=30 \
       binsz=0.05 \
       coordsys="CEL" \
       xref=86.63 \
       yref=22.01 \
       proj="CAR"
$ECHO -n "."
if [ -s "cntmap2.fits" ]
then
  $ECHO -n "."
else
  $ECHO " cntmap2.fits file is not found"
  exit 1
fi
$ECHO " ok"

#
# Test ctexpcube
# ==============
$ECHO -n "Test ctexpcube: "
#
# Run 1
$ctexpcube inobs="data/crab_events.fits" \
           incube="NONE" \
           outcube="expcube1.fits" \
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
$ECHO " expcube1.fits file is not found"
exit 1
fi
#
# Run 2
$ctexpcube inobs="events.fits" \
           incube="cntmap2.fits" \
           outcube="expcube2.fits" \
           caldb="irf" \
           irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "expcube2.fits" ]
then
$ECHO -n "."
else
$ECHO " expcube2.fits file is not found"
exit 1
fi
$ECHO " ok"


#
# Test ctpsfcube
# ==============
$ECHO -n "Test ctpsfcube: "
#
# Run 1
$ctpsfcube inobs="data/crab_events.fits" \
           incube="NONE" \
           outcube="psfcube1.fits" \
           caldb="irf" \
           irf="cta_dummy_irf" \
           emin=0.1 \
           emax=100.0 \
           enumbins=20 \
           nxpix=10 \
           nypix=10 \
           binsz=0.4 \
           coordsys="CEL" \
           xref=83.63 \
           yref=22.01 \
           proj="CAR" \
           amax=0.3 \
           anumbins=10
$ECHO -n "."
if [ -s "psfcube1.fits" ]
then
$ECHO -n "."
else
$ECHO " psfcube1.fits file is not found"
exit 1
fi
#
# Run 2
$ctpsfcube inobs="events.fits" \
           incube="cntmap2.fits" \
           outcube="psfcube2.fits" \
           caldb="irf" \
           irf="cta_dummy_irf" \
           amax=0.3 \
           anumbins=2
$ECHO -n "."
if [ -s "psfcube2.fits" ]
then
$ECHO -n "."
else
$ECHO " psfcube2.fits file is not found"
exit 1
fi
$ECHO " ok"


#
# Test ctbkgcube
# ==============
$ECHO -n "Test ctbkgcube: "
#
# Run 1
$ctbkgcube inobs="data/crab_events.fits" \
           inmodel="data/crab.xml" \
           incube="NONE" \
           outcube="bkgcube1.fits" \
           outmodel="bkgcube1.xml" \
           caldb="irf" \
           irf="cta_dummy_irf" \
           emin=0.1 \
           emax=100.0 \
           enumbins=20 \
           nxpix=10 \
           nypix=10 \
           binsz=0.4 \
           coordsys="CEL" \
           xref=83.63 \
           yref=22.01 \
           proj="CAR"
$ECHO -n "."
if [ -s "bkgcube1.fits" ]
then
$ECHO -n "."
else
$ECHO " bkgcube1.fits file is not found"
exit 1
fi
if [ -s "bkgcube1.xml" ]
then
$ECHO -n "."
else
$ECHO " bkgcube1.xml file is not found"
exit 1
fi
#
# Run 2
$ctbkgcube inobs="events.fits" \
           inmodel="data/crab.xml" \
           incube="cntmap2.fits" \
           caldb="irf" \
           irf="cta_dummy_irf" \
           outcube="bkgcube2.fits" \
           outmodel="bkgcube2.xml"
$ECHO -n "."
if [ -s "bkgcube2.fits" ]
then
$ECHO -n "."
else
$ECHO " bkgcube2.fits file is not found"
exit 1
fi
if [ -s "bkgcube2.fits" ]
then
$ECHO -n "."
else
$ECHO " bkgcube2.xml file is not found"
exit 1
fi
$ECHO " ok"



#
# Test ctmodel
# ============
$ECHO -n "Test ctmodel: "
#
# Run 1
$ctmodel inobs="cntmap1.fits" \
         incube="cntmap1.fits" \
         outcube="modmap1.fits" \
         expcube="NONE" \
         psfcube="NONE" \
         bkgcube="NONE" \
         caldb="irf" \
         irf="cta_dummy_irf" \
         inmodel="data/crab.xml"
$ECHO -n "."
if [ -s "modmap1.fits" ]
then
  $ECHO -n "."
else
  $ECHO " modmap1.fits file is not found"
  exit 1
fi

#
# Run 2
$ctmodel inobs="cntmap2.fits" \
         incube="cntmap2.fits" \
         outcube="modmap2.fits" \
         expcube="expcube2.fits" \
         psfcube="psfcube2.fits" \
         bkgcube="bkgcube2.fits" \
         caldb="irf" \
         irf="cta_dummy_irf" \
         inmodel="bkgcube2.xml"
$ECHO -n "."
if [ -s "modmap2.fits" ]
then
$ECHO -n "."
else
$ECHO " modmap2.fits file is not found"
exit 1
fi

#
# Run 3
$ctmodel inobs="NONE" \
         incube="NONE" \
         outcube="modmap3.fits" \
         expcube="NONE" \
         psfcube="NONE" \
         bkgcube="NONE" \
         caldb="irf" \
         irf="cta_dummy_irf" \
         inmodel="data/crab.xml" \
         ra=83.63 \
         dec=22.01 \
         rad=10.0 \
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
  $ECHO " modmap3.fits file is not found"
  exit 1
fi
$ECHO " ok"


#
# Test ctselect
# =============
$ECHO -n "Test ctselect: "
$ctselect inobs="events.fits" \
          outobs="selected_events.fits" \
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
  $ECHO " selected_events.fits file is not found"
  exit 1
fi
$ECHO " ok"


#
# Test ctlike
# ===========
$ECHO -n "Test ctlike: "
#
# Run 1
$ctlike inobs="cntmap1.fits" \
        inmodel="data/crab.xml" \
        outmodel="results_binned.xml" \
        expcube="NONE" \
        psfcube="NONE" \
        bkgcube="NONE" \
        caldb="irf" \
        irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "results_binned.xml" ]
then
  $ECHO -n "."
else
  $ECHO " results_binned.xml file is not found"
  exit 1
fi
#
# Run 2
$ctlike inobs="selected_events.fits" \
        inmodel="data/crab.xml" \
        outmodel="results_unbinned.xml" \
        caldb="irf" \
        irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "results_unbinned.xml" ]
then
  $ECHO -n "."
else
  $ECHO " results_unbinned.xml file is not found"
  exit 1
fi

#
# Run 3
$ctlike inobs="cntmap2.fits" \
        inmodel="bkgcube2.xml" \
        outmodel="results_binned_cube_background.xml" \
        expcube="expcube2.fits" \
        psfcube="psfcube2.fits" \
        bkgcube="bkgcube2.fits" \
        caldb="irf" \
        irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "results_binned_cube_background.xml" ]
then
  $ECHO -n "."
else
  $ECHO " results_binned_cube_background.xml file is not found"
  exit 1
fi
$ECHO " ok"


#
# Test cttsmap
# ============
$ECHO -n "Test cttsmap: "
$cttsmap inobs="selected_events.fits" \
         inmodel="data/crab.xml" \
         srcname="Crab" \
         caldb="irf" \
         irf="cta_dummy_irf" \
         outmap="tsmap.fits" \
         nxpix=5 \
         nypix=5 \
         binsz=0.05 \
         coordsys="CEL" \
         xref=83.63 \
         yref=22.01 \
         proj="CAR"
$ECHO -n "."
if [ -s "tsmap.fits" ]
then
  $ECHO -n "."
else
  $ECHO " tsmap.fits file is not found"
  exit 1
fi
$ECHO " ok"




#
# Test ctcubemask
# ===============
$ECHO -n "Test ctcubemask: "
$ctcubemask inobs="data/crab_cntmap.fits" \
            regfile="data/exclusion.reg" \
            outcube="filtered_cube.fits" \
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
  $ECHO " filtered_cube.fits file is not found"
  exit 1
fi
$ECHO " ok"


#
# Test ctbutterfly
# ===============
$ECHO -n "Test ctbutterfly: "
$ctbutterfly inobs="data/crab_events.fits.gz" \
             inmodel="data/crab.xml" \
             srcname="Crab" \
             outfile="butterfly.txt" \
             caldb="irf" \
             irf="cta_dummy_irf" \
             emin=0.1 \
             emax=100.0
$ECHO -n "."
if [ -s "butterfly.txt" ]
then
    $ECHO -n "."
else
    $ECHO " butterfly.txt file is not found"
    exit 1
fi
$ECHO " ok"


#
# Test ctulimit
# ===============
$ECHO -n "Test ctulimit: "
$ctulimit inobs="data/crab_events.fits.gz" \
          inmodel="data/crab.xml" \
          srcname="Crab" \
          caldb="irf" \
          irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "ctulimit.log" ]
then
    $ECHO -n "."
else
    $ECHO " ctulimit.log file is not found"
    exit 1
fi
$ECHO " ok"


#
# Test cterror
# ============
$ECHO -n "Test cterror: "
$cterror inobs="data/crab_events.fits.gz" \
         inmodel="data/crab.xml" \
         outmodel="results_error.xml" \
         srcname="Crab" \
         caldb="irf" \
         irf="cta_dummy_irf"
$ECHO -n "."
if [ -s "cterror.log" ]
then
    $ECHO -n "."
else
    $ECHO " cterror.log file is not found"
    exit 1
fi
$ECHO " ok"

