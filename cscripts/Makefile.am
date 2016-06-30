## Process this file with automake to produce Makefile.in

# Scripts to be distributed in $(pyexecdir)/cscripts
pythonscripts = __init__.py \
                cscaldb.py \
                csinfo.py \
                cslightcrv.py \
                csmodelinfo.py \
                csmodelmerge.py \
                csobs2caldb.py \
                csobsdef.py \
                csobsinfo.py \
                cspull.py \
                csresmap.py \
                csroot2caldb.py  \
                cssens.py \
                csspec.py \
                cstsdist.py \
                cstsmapmerge.py \
                csworkflow.py \
                csfindobs.py \
                csiactcopy.py \
                csiactdata.py \
                csiactobs.py \
                cstssplit.py \
                obsutils.py

# Parfiles to be distributed in $(prefix)/syspfiles
parfiles = cscaldb.par \
           cslightcrv.par \
           csmodelinfo.par \
           csmodelmerge.par \
           csobs2caldb.par \
           csobsdef.par \
           csobsinfo.par \
           cspull.par \
           csresmap.par \
           csroot2caldb.par \
           cssens.par \
           csspec.par \
           cstsdist.par \
           cstsmapmerge.par \
           csworkflow.par \
           csfindobs.par \
           csiactcopy.par \
           csiactdata.par \
           cstssplit.par \
           csiactobs.par

# Scripts for which symbolic links will be created in bindir
linkscripts = cscaldb \
              csinfo \
              cslightcrv \
              csmodelinfo \
              csmodelmerge \
              csobs2caldb \
              csobsdef \
              csobsinfo \
              cspull \
              csresmap \
              csroot2caldb \
              cssens \
              csspec \
              cstsdist \
              cstsmapmerge \
              csworkflow \
              csfindobs \
              csiactcopy \
              csiactdata \
              cstssplit \
              csiactobs

# Test scripts to be distributed in $(pyexecdir)/cscripts/tests
testscripts = tests/__init__.py \
              ../test/test_python_cscripts.py \
              ../test/test_cscaldb.py \
              ../test/test_csmodelinfo.py \
              ../test/test_csmodelmerge.py \
              ../test/test_csobs2caldb.py \
              ../test/test_csobsdef.py \
              ../test/test_csobsinfo.py \
              ../test/test_cslightcrv.py \
              ../test/test_cspull.py \
              ../test/test_csresmap.py \
              ../test/test_csroot2caldb.py \
              ../test/test_cssens.py \
              ../test/test_csspec.py \
              ../test/test_cstsdist.py \
              ../test/test_cstsmapmerge.py \
              ../test/test_csworkflow.py \
              ../test/test_csfindobs.py \
              ../test/test_csiactcopy.py \
              ../test/test_csiactdata.py \
              ../test/test_cstssplit.py \
              ../test/test_csiactobs.py

# Test data to be distributed in $(pyexecdir)/cscripts/tests/data
testdata    = ../test/data/crab.xml \
              ../test/data/obs_unbinned.xml \
              ../test/data/obs_stacked.xml \
              ../test/data/crab_events.fits.gz \
              ../test/data/crab_cntmap.fits.gz \
              ../test/data/crab_cntmap_small.fits.gz \
              ../test/data/crab_expcube.fits.gz \
              ../test/data/crab_psfcube.fits.gz \
              ../test/data/crab_edispcube.fits.gz \
              ../test/data/crab_bkgcube.fits.gz \
              ../test/data/tsmap_0.fits \
              ../test/data/tsmap_1.fits \
              ../test/data/exclusion.reg \
              ../test/data/pntdef.dat \
              ../test/data/models.txt \
              ../test/data/tsmaps.txt \
              ../test/data/model_cube_background.xml \
              ../test/data/model_cube_background2.xml \
              ../test/data/crab_bkgcube.xml \
              ../test/data/irf.root \
              ../test/data/workflow.xml

# Test data to be distributed in $(pyexecdir)/cscripts/tests/iactdata
testiactdata = ../test/iactdata/hdu-index.fits.gz \
               ../test/iactdata/master.json \
               ../test/iactdata/obs-index.fits.gz

# Test IRFs to be distributed in $(pyexecdir)/cscripts/tests/irf
testirf      = ../test/irf/irf_file.fits.gz \
               ../test/irf/cta_dummy_irf.dat

# Distribute scripts and parfiles
EXTRA_DIST         = $(binscripts) $(pythonscripts)
binscriptsdir      = $(bindir)
binscripts_DATA    = $(binscripts)
pythonscriptsdir   = $(pyexecdir)/cscripts
pythonscripts_DATA = $(pythonscripts)
pfilesdir          = $(prefix)/syspfiles
pfiles_DATA        = $(parfiles)
testscriptsdir     = $(pyexecdir)/cscripts/tests
testscripts_DATA   = $(testscripts)
testdatadir        = $(pyexecdir)/cscripts/tests/data
testdata_DATA      = $(testdata)
testiactdatadir    = $(pyexecdir)/cscripts/tests/iactdata
testiactdata_DATA  = $(testiactdata)
testirfdir         = $(pyexecdir)/cscripts/tests/irf
testirf_DATA       = $(testirf)

# Add hook rule to create symbolic links for cscripts
install-data-hook:
	@binlist='$(binscripts)';\
	for p in $$binlist; do \
	  chmod +x "$(DESTDIR)$(bindir)/$$p"; \
	done
	@linklist='$(linkscripts)';\
	for p in $$linklist; do \
      rm -rf "$(DESTDIR)$(bindir)/$$p"; \
      ln -s "$(DESTDIR)$(pyexecdir)/cscripts/$$p.py" "$(DESTDIR)$(bindir)/$$p"; \
	  chmod +x "$(DESTDIR)$(pyexecdir)/cscripts/$$p.py"; \
	done
    
