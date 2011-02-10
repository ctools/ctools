#! /bin/bash
#
# ctatools test script
# ====================================================================

#
# Print Header
#
echo "*****************************"
echo "* Test ctatools executables *"
echo "*****************************"


#
# Remove any existing result files
# ================================
rm -rf *.fits *.log *.xml


#
# Test ctobssim
# ============
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
