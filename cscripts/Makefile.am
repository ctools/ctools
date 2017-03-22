#############################################################################
# Makefile for cscripts scripts                                             #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2017 Juergen Knoedlseder                               #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
#############################################################################

# Scripts to be distributed in $(pyexecdir)/cscripts
pythonscripts = $(top_builddir)/cscripts/__init__.py \
                $(srcdir)/cscaldb.py \
                $(srcdir)/csebins.py \
                $(srcdir)/csinfo.py \
                $(srcdir)/cslightcrv.py \
                $(srcdir)/csmodelinfo.py \
                $(srcdir)/csmodelmerge.py \
                $(srcdir)/csmodelselect.py \
                $(srcdir)/csobs2caldb.py \
                $(srcdir)/csobsdef.py \
                $(srcdir)/csobsinfo.py \
                $(srcdir)/csobsselect.py \
                $(srcdir)/cspull.py \
                $(srcdir)/csresmap.py \
                $(srcdir)/csroot2caldb.py  \
                $(srcdir)/cssens.py \
                $(srcdir)/csspec.py \
                $(srcdir)/cssrcdetect.py \
                $(srcdir)/cstsdist.py \
                $(srcdir)/cstsmapmerge.py \
                $(srcdir)/cstsmapsplit.py \
                $(srcdir)/csviscube.py \
                $(srcdir)/csworkflow.py \
                $(srcdir)/csfindobs.py \
                $(srcdir)/csiactcopy.py \
                $(srcdir)/csiactdata.py \
                $(srcdir)/csiactobs.py \
                $(srcdir)/obsutils.py \
                $(srcdir)/modutils.py \
                $(srcdir)/calutils.py \
                $(srcdir)/ioutils.py

# Parfiles to be distributed in $(prefix)/syspfiles
parfiles = $(srcdir)/cscaldb.par \
           $(srcdir)/csebins.par \
           $(srcdir)/cslightcrv.par \
           $(srcdir)/csmodelinfo.par \
           $(srcdir)/csmodelmerge.par \
           $(srcdir)/csmodelselect.par \
           $(srcdir)/csobs2caldb.par \
           $(srcdir)/csobsdef.par \
           $(srcdir)/csobsinfo.par \
           $(srcdir)/csobsselect.par \
           $(srcdir)/cspull.par \
           $(srcdir)/csresmap.par \
           $(srcdir)/csroot2caldb.par \
           $(srcdir)/cssens.par \
           $(srcdir)/csspec.par \
           $(srcdir)/cssrcdetect.par \
           $(srcdir)/cstsdist.par \
           $(srcdir)/cstsmapmerge.par \
           $(srcdir)/cstsmapsplit.par \
           $(srcdir)/csviscube.par \
           $(srcdir)/csworkflow.par \
           $(srcdir)/csfindobs.par \
           $(srcdir)/csiactcopy.par \
           $(srcdir)/csiactdata.par \
           $(srcdir)/csiactobs.par

# Scripts for which symbolic links will be created in bindir
linkscripts = cscaldb \
              csebins \
              csinfo \
              cslightcrv \
              csmodelinfo \
              csmodelmerge \
              csmodelselect \
              csobs2caldb \
              csobsdef \
              csobsinfo \
              csobsselect \
              cspull \
              csresmap \
              csroot2caldb \
              cssens \
              csspec \
              cssrcdetect \
              cstsdist \
              cstsmapmerge \
              cstsmapsplit \
              csviscube \
              csworkflow \
              csfindobs \
              csiactcopy \
              csiactdata \
              csiactobs

# Test scripts to be distributed in $(pyexecdir)/cscripts/tests
testscripts = $(srcdir)/tests/__init__.py \
              $(top_srcdir)/test/testing.py \
              $(top_srcdir)/test/test_python_cscripts.py \
              $(top_srcdir)/test/test_cscript.py \
              $(top_srcdir)/test/test_cscaldb.py \
              $(top_srcdir)/test/test_csebins.py \
              $(top_srcdir)/test/test_csfindobs.py \
              $(top_srcdir)/test/test_csinfo.py \
              $(top_srcdir)/test/test_cslightcrv.py \
              $(top_srcdir)/test/test_csmodelinfo.py \
              $(top_srcdir)/test/test_csmodelmerge.py \
              $(top_srcdir)/test/test_csmodelselect.py \
              $(top_srcdir)/test/test_csobs2caldb.py \
              $(top_srcdir)/test/test_csobsdef.py \
              $(top_srcdir)/test/test_csobsinfo.py \
              $(top_srcdir)/test/test_csobsselect.py \
              $(top_srcdir)/test/test_cspull.py \
              $(top_srcdir)/test/test_csresmap.py \
              $(top_srcdir)/test/test_csroot2caldb.py \
              $(top_srcdir)/test/test_cssens.py \
              $(top_srcdir)/test/test_csspec.py \
              $(top_srcdir)/test/test_cssrcdetect.py \
              $(top_srcdir)/test/test_cstsdist.py \
              $(top_srcdir)/test/test_cstsmapmerge.py \
              $(top_srcdir)/test/test_cstsmapsplit.py \
              $(top_srcdir)/test/test_csviscube.py \
              $(top_srcdir)/test/test_csworkflow.py \
              $(top_srcdir)/test/test_csiactcopy.py \
              $(top_srcdir)/test/test_csiactdata.py \
              $(top_srcdir)/test/test_csiactobs.py

# Test data to be distributed in $(pyexecdir)/cscripts/tests/data
testdata = $(top_srcdir)/test/data/crab.xml \
           $(top_srcdir)/test/data/obs_unbinned.xml \
           $(top_srcdir)/test/data/obs_stacked.xml \
           $(top_srcdir)/test/data/obs_stacked_two.xml \
           $(top_srcdir)/test/data/crab_events.fits.gz \
           $(top_srcdir)/test/data/crab_cntmap.fits.gz \
           $(top_srcdir)/test/data/crab_cntmap_small.fits.gz \
           $(top_srcdir)/test/data/crab_expcube.fits.gz \
           $(top_srcdir)/test/data/crab_psfcube.fits.gz \
           $(top_srcdir)/test/data/crab_edispcube.fits.gz \
           $(top_srcdir)/test/data/crab_bkgcube.fits.gz \
           $(top_srcdir)/test/data/tsmap_0.fits \
           $(top_srcdir)/test/data/tsmap_1.fits \
           $(top_srcdir)/test/data/exclusion.reg \
           $(top_srcdir)/test/data/pntdef_minimal.dat \
           $(top_srcdir)/test/data/pntdef_maximal.dat \
           $(top_srcdir)/test/data/models.txt \
           $(top_srcdir)/test/data/tsmaps.txt \
           $(top_srcdir)/test/data/model_cube_background.xml \
           $(top_srcdir)/test/data/model_cube_background2.xml \
           $(top_srcdir)/test/data/crab_bkgcube.xml \
           $(top_srcdir)/test/data/irf.root \
           $(top_srcdir)/test/data/irf_file.fits.gz \
           $(top_srcdir)/test/data/workflow.xml \
           $(top_srcdir)/test/data/butterfly.txt \
           $(top_srcdir)/test/data/pha.fits \
           $(top_srcdir)/test/data/pha_channels.fits \
           $(top_srcdir)/test/data/pull.dat \
           $(top_srcdir)/test/data/sensitivity.dat \
           $(top_srcdir)/test/data/sensitivity_integral.dat \
           $(top_srcdir)/test/data/skymap.fits \
           $(top_srcdir)/test/data/spectrum.fits \
           $(top_srcdir)/test/data/lightcurve.fits \
           $(top_srcdir)/test/data/ts.dat \
           $(top_srcdir)/test/data/models_spatial.xml \
           $(top_srcdir)/test/data/models_spectral.xml \
           $(top_srcdir)/test/data/model_filefunction.txt \
           $(top_srcdir)/test/data/model_map.fits.gz \
           $(top_srcdir)/test/data/model_mapcube.fits.gz \
           $(top_srcdir)/test/data/iact_runlist.dat \
           $(top_srcdir)/test/data/iact_runlist_1.dat \
           $(top_srcdir)/test/data/iact_runlist_2.dat

# Test data to be distributed in $(pyexecdir)/cscripts/tests/data/iactdata
testiactdata = $(top_srcdir)/test/data/iactdata/hdu-index.fits.gz \
               $(top_srcdir)/test/data/iactdata/master.json \
               $(top_srcdir)/test/data/iactdata/obs-index.fits.gz

# Test data to be distributed in $(pyexecdir)/cscripts/tests/data/iactdata/events
testiactevents = $(top_srcdir)/test/data/iactdata/events/irf_file.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_0.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_1.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_2.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_3.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_4.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_5.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_6.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_7.fits.gz \
                 $(top_srcdir)/test/data/iactdata/events/events_8.fits.gz

# Install scripts, parfiles, test scripts and test data. Note that we don't
# need to prepend $(DESTDIR) here (not fully clear why, but that's the way
# it is)
binscriptsdir       = $(bindir)
binscripts_DATA     = $(binscripts)
pythonscriptsdir    = $(pyexecdir)/cscripts
pythonscripts_DATA  = $(pythonscripts)
pfilesdir           = $(prefix)/syspfiles
pfiles_DATA         = $(parfiles)
testscriptsdir      = $(pyexecdir)/cscripts/tests
testscripts_DATA    = $(testscripts)
testdatadir         = $(pyexecdir)/cscripts/tests/data
testdata_DATA       = $(testdata)
testiactdatadir     = $(pyexecdir)/cscripts/tests/data/iactdata
testiactdata_DATA   = $(testiactdata)
testiacteventsdir   = $(pyexecdir)/cscripts/tests/data/iactdata/events
testiactevents_DATA = $(testiactevents)

# Files to be distributed in the tar file
EXTRA_DIST = $(binscripts) $(pythonscripts) $(parfiles) \
             $(testscripts) $(testdata) $(testiactdata) \
             $(testiactevents)

# Add install hook rule to create symbolic links for cscripts
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

# Uninstall hook rule to remove the symbolic links for cscripts
uninstall-hook:
	@linklist='$(linkscripts)';\
	for p in $$linklist; do \
		rm -rf "$(DESTDIR)$(bindir)/$$p"; \
	done

# Files to be removed in build directory
clean-local:
	@if ! [ '@top_srcdir@' -ef '@top_builddir@' ]; then \
		rm -rf $(top_builddir)/cscripts/cs*.py; \
		rm -rf $(top_builddir)/cscripts/obsutils.py; \
		rm -rf $(top_builddir)/cscripts/modutils.py; \
		rm -rf $(top_builddir)/cscripts/calutils.py; \
		rm -rf $(top_builddir)/cscripts/ioutils.py; \
		rm -rf $(top_builddir)/cscripts/*.par; \
	fi
	rm -rf $(top_builddir)/cscripts/*.pyc
	rm -rf $(top_builddir)/cscripts/__pycache__

# Uninstall the cscripts Python module
uninstall-local:
	rm -rf $(DESTDIR)$(pyexecdir)/cscripts/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/cscripts/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/cscripts/__pycache__
	rm -rf $(DESTDIR)$(pyexecdir)/cscripts/tests/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/cscripts/tests/__pycache__
