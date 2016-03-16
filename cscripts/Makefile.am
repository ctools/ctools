## Process this file with automake to produce Makefile.in

# Scripts to be distributed in $(pyexecdir)/cscripts
pythonscripts = __init__.py \
                cspull.py cssens.py cstsdist.py cscaldb.py csresmap.py \
                csobsdef.py csspec.py cslightcrv.py csobsinfo.py \
                csinfo.py csmodelinfo.py csmodelmerge.py cstsmapmerge.py \
                csobs2caldb.py \
                csfindobs.py csiactcopy.py csiactdata.py csiactobs.py \
                obsutils.py

# Parfiles to be distributed in $(prefix)/syspfiles
parfiles = cspull.par cssens.par cstsdist.par cscaldb.par csresmap.par \
           csobsdef.par csspec.par cslightcrv.par csobsinfo.par \
           csmodelinfo.par csmodelmerge.par cstsmapmerge.par \
           csobs2caldb.par \
           csfindobs.par csiactcopy.par csiactdata.par csiactobs.par

# Scripts for which symbolic links will be created in bindir
linkscripts = cspull cssens cstsdist cscaldb csresmap \
              csobsdef csspec cslightcrv csobsinfo \
              csinfo csmodelinfo csmodelmerge cstsmapmerge \
              csobs2caldb \
              csfindobs csiactcopy csiactdata csiactobs

# Test scripts to be distributed in $(pyexecdir)/cscripts/tests
testscripts = tests/__init__.py \
              ../test/test_python_cscripts.py \
              ../test/test_cscaldb.py \
              ../test/test_csobs2caldb.py \
              ../test/test_cslightcrv.py

# Test data to be distributed in $(pyexecdir)/cscripts/tests/data
testdata    = ../test/data/crab.xml \
              ../test/data/obs_unbinned.xml \
              ../test/data/crab_events.fits.gz \
              ../test/data/crab_cntmap.fits.gz \
              ../test/data/exclusion.reg

# Test IRFs to be distributed in $(pyexecdir)/cscripts/tests/irf
#testirf     = 

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
#testirfdir         = $(pyexecdir)/cscripts/tests/irf
#testirf_DATA       = $(testirf)

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
    
