#############################################################################
# Makefile for cscripts scripts                                             #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2022 Juergen Knoedlseder                               #
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
                $(srcdir)/csadd2caldb.py \
                $(srcdir)/csbkgmodel.py \
                $(srcdir)/cscaldb.py \
                $(srcdir)/csebins.py \
                $(srcdir)/csfootprint.py \
                $(srcdir)/csinfo.py \
                $(srcdir)/cslightcrv.py \
                $(srcdir)/csmodelinfo.py \
                $(srcdir)/csmodelmerge.py \
                $(srcdir)/csmodelselect.py \
                $(srcdir)/csmodelsois.py \
                $(srcdir)/csobs2caldb.py \
                $(srcdir)/csobsdef.py \
                $(srcdir)/csobsinfo.py \
                $(srcdir)/csobsselect.py \
                $(srcdir)/csphagen.py \
                $(srcdir)/csphasecrv.py \
                $(srcdir)/cspull.py \
                $(srcdir)/csresmap.py \
                $(srcdir)/csresspec.py \
                $(srcdir)/csroot2caldb.py  \
                $(srcdir)/cssens.py \
                $(srcdir)/csspec.py \
                $(srcdir)/cssrcdetect.py \
                $(srcdir)/cstsdist.py \
                $(srcdir)/cstsmapmerge.py \
                $(srcdir)/cstsmapsplit.py \
                $(srcdir)/csviscube.py \
                $(srcdir)/csscs.py \
                $(srcdir)/csworkflow.py \
                $(srcdir)/csfindobs.py \
                $(srcdir)/csiactcopy.py \
                $(srcdir)/csiactdata.py \
                $(srcdir)/csiactobs.py \
                $(srcdir)/obsutils.py \
                $(srcdir)/modutils.py \
                $(srcdir)/calutils.py \
                $(srcdir)/ioutils.py \
                $(srcdir)/mputils.py 	

# Parfiles to be distributed in $(prefix)/syspfiles
parfiles = $(srcdir)/csadd2caldb.par \
           $(srcdir)/csbkgmodel.par \
           $(srcdir)/cscaldb.par \
           $(srcdir)/csebins.par \
           $(srcdir)/cslightcrv.par \
           $(srcdir)/csmodelinfo.par \
           $(srcdir)/csmodelmerge.par \
           $(srcdir)/csmodelselect.par \
           $(srcdir)/csmodelsois.par \
           $(srcdir)/csobs2caldb.par \
           $(srcdir)/csobsdef.par \
           $(srcdir)/csobsinfo.par \
           $(srcdir)/csobsselect.par \
           $(srcdir)/csphagen.par \
           $(srcdir)/csphasecrv.par \
           $(srcdir)/cspull.par \
           $(srcdir)/csresmap.par \
           $(srcdir)/csresspec.par \
           $(srcdir)/csroot2caldb.par \
           $(srcdir)/cssens.par \
           $(srcdir)/csspec.par \
           $(srcdir)/cssrcdetect.par \
           $(srcdir)/cstsdist.par \
           $(srcdir)/cstsmapmerge.par \
           $(srcdir)/cstsmapsplit.par \
           $(srcdir)/csviscube.par \
           $(srcdir)/csscs.par \
           $(srcdir)/csfootprint.par \
           $(srcdir)/csworkflow.par \
           $(srcdir)/csfindobs.par \
           $(srcdir)/csiactcopy.par \
           $(srcdir)/csiactdata.par \
           $(srcdir)/csiactobs.par

# Scripts for which symbolic links will be created in bindir
linkscripts = csadd2caldb \
              csbkgmodel \
              cscaldb \
              csebins \
              csinfo \
              cslightcrv \
              csmodelinfo \
              csmodelmerge \
              csmodelselect \
              csmodelsois \
              csobs2caldb \
              csobsdef \
              csobsinfo \
              csobsselect \
              csphagen \
              csphasecrv \
              cspull \
              csresmap \
              csresspec \
              csroot2caldb \
              cssens \
              csspec \
              cssrcdetect \
              cstsdist \
              cstsmapmerge \
              cstsmapsplit \
              csviscube \
              csscs \
              csfootprint \
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
              $(top_srcdir)/test/test_csadd2caldb.py \
              $(top_srcdir)/test/test_csbkgmodel.py \
              $(top_srcdir)/test/test_csebins.py \
              $(top_srcdir)/test/test_csfindobs.py \
              $(top_srcdir)/test/test_csinfo.py \
              $(top_srcdir)/test/test_cslightcrv.py \
              $(top_srcdir)/test/test_csmodelinfo.py \
              $(top_srcdir)/test/test_csmodelmerge.py \
              $(top_srcdir)/test/test_csmodelselect.py \
              $(top_srcdir)/test/test_csmodelsois.py \
              $(top_srcdir)/test/test_csobs2caldb.py \
              $(top_srcdir)/test/test_csobsdef.py \
              $(top_srcdir)/test/test_csobsinfo.py \
              $(top_srcdir)/test/test_csobsselect.py \
              $(top_srcdir)/test/test_csphagen.py \
              $(top_srcdir)/test/test_csphasecrv.py \
              $(top_srcdir)/test/test_cspull.py \
              $(top_srcdir)/test/test_csresmap.py \
              $(top_srcdir)/test/test_csresspec.py \
              $(top_srcdir)/test/test_csroot2caldb.py \
              $(top_srcdir)/test/test_cssens.py \
              $(top_srcdir)/test/test_csspec.py \
              $(top_srcdir)/test/test_cssrcdetect.py \
              $(top_srcdir)/test/test_cstsdist.py \
              $(top_srcdir)/test/test_cstsmapmerge.py \
              $(top_srcdir)/test/test_cstsmapsplit.py \
              $(top_srcdir)/test/test_csviscube.py \
              $(top_srcdir)/test/test_csscs.py \
              $(top_srcdir)/test/test_csfootprint.py \
              $(top_srcdir)/test/test_csworkflow.py \
              $(top_srcdir)/test/test_obsutils.py \
              $(top_srcdir)/test/test_mputils.py \
              $(top_srcdir)/test/test_csiactcopy.py \
              $(top_srcdir)/test/test_csiactdata.py \
              $(top_srcdir)/test/test_csiactobs.py

# Test data to be distributed in $(pyexecdir)/cscripts/tests/data
testdata = $(top_srcdir)/test/data/crab.xml \
           $(top_srcdir)/test/data/obs_unbinned.xml \
           $(top_srcdir)/test/data/obs_unbinned_two.xml \
           $(top_srcdir)/test/data/obs_stacked.xml \
           $(top_srcdir)/test/data/obs_stacked_two.xml \
           $(top_srcdir)/test/data/crab_events.fits \
           $(top_srcdir)/test/data/crab_cntmap.fits \
           $(top_srcdir)/test/data/crab_expcube.fits.gz \
           $(top_srcdir)/test/data/crab_psfcube.fits.gz \
           $(top_srcdir)/test/data/crab_edispcube.fits.gz \
           $(top_srcdir)/test/data/crab_bkgcube.fits.gz \
           $(top_srcdir)/test/data/crab_onoff.xml \
           $(top_srcdir)/test/data/crab_offaxis1.fits.gz \
           $(top_srcdir)/test/data/crab_offaxis2.fits.gz \
           $(top_srcdir)/test/data/crab_exclusion.fits.gz \
           $(top_srcdir)/test/data/background_events.fits.gz \
           $(top_srcdir)/test/data/crab_src_reg.reg \
           $(top_srcdir)/test/data/crab_bkg_reg.reg \
           $(top_srcdir)/test/data/tsmap_0.fits \
           $(top_srcdir)/test/data/tsmap_1.fits \
           $(top_srcdir)/test/data/exclusion.reg \
           $(top_srcdir)/test/data/pntdef_minimal.dat \
           $(top_srcdir)/test/data/pntdef_maximal.dat \
           $(top_srcdir)/test/data/models.txt \
           $(top_srcdir)/test/data/tsmaps.txt \
           $(top_srcdir)/test/data/model_cube_background1.xml \
           $(top_srcdir)/test/data/model_cube_background2.xml \
           $(top_srcdir)/test/data/crab_bkgcube.xml \
           $(top_srcdir)/test/data/irf.root \
           $(top_srcdir)/test/data/irf_file.fits.gz \
           $(top_srcdir)/test/data/workflow.xml \
           $(top_srcdir)/test/data/butterfly.txt \
           $(top_srcdir)/test/data/butterfly.fits \
           $(top_srcdir)/test/data/pha.fits \
           $(top_srcdir)/test/data/pha_channels.fits \
           $(top_srcdir)/test/data/pull.dat \
           $(top_srcdir)/test/data/pull.fits \
           $(top_srcdir)/test/data/sensitivity.dat \
           $(top_srcdir)/test/data/sensitivity_integral.dat \
           $(top_srcdir)/test/data/sensitivity.fits \
           $(top_srcdir)/test/data/skymap.fits \
           $(top_srcdir)/test/data/spectrum.fits \
           $(top_srcdir)/test/data/spectrum_castro.fits \
           $(top_srcdir)/test/data/statistics.xml \
           $(top_srcdir)/test/data/residual.fits \
           $(top_srcdir)/test/data/residual_onoff.fits \
           $(top_srcdir)/test/data/lightcurve.fits \
           $(top_srcdir)/test/data/ts.dat \
           $(top_srcdir)/test/data/ts.fits \
           $(top_srcdir)/test/data/models_spatial.xml \
           $(top_srcdir)/test/data/models_spectral.xml \
           $(top_srcdir)/test/data/model_filefunction.txt \
           $(top_srcdir)/test/data/model_map.fits.gz \
           $(top_srcdir)/test/data/model_mapcube.fits.gz \
           $(top_srcdir)/test/data/iact_runlist.dat \
           $(top_srcdir)/test/data/iact_runlist_1.dat \
           $(top_srcdir)/test/data/iact_runlist_2.dat \
           $(top_srcdir)/test/data/model_temporal_phasecurve.xml \
           $(top_srcdir)/test/data/model_temporal_phasecurve.fits \
           $(top_srcdir)/test/data/phased_events.fits.gz \
           $(top_srcdir)/test/data/onoff_arf.fits.gz \
           $(top_srcdir)/test/data/onoff_arf1.fits.gz \
           $(top_srcdir)/test/data/onoff_arf2.fits.gz \
           $(top_srcdir)/test/data/onoff_model.xml \
           $(top_srcdir)/test/data/onoff_obs.xml \
           $(top_srcdir)/test/data/onoff_obs1.xml \
           $(top_srcdir)/test/data/onoff_obs2.xml \
           $(top_srcdir)/test/data/onoff_obs_stacked.xml \
           $(top_srcdir)/test/data/onoff_pha_off.fits.gz \
           $(top_srcdir)/test/data/onoff_pha_off1.fits.gz \
           $(top_srcdir)/test/data/onoff_pha_off2.fits.gz \
           $(top_srcdir)/test/data/onoff_pha_on.fits.gz \
           $(top_srcdir)/test/data/onoff_pha_on1.fits.gz \
           $(top_srcdir)/test/data/onoff_pha_on2.fits.gz \
           $(top_srcdir)/test/data/onoff_rmf.fits.gz \
           $(top_srcdir)/test/data/onoff_rmf1.fits.gz \
           $(top_srcdir)/test/data/onoff_rmf2.fits.gz \
           $(top_srcdir)/test/data/model_crab_radialacceptance.xml \
           $(top_srcdir)/test/data/sigdist_map.fits.gz \
           $(top_srcdir)/test/data/variability.fits.gz \
           $(top_srcdir)/test/data/viscube.fits.gz \
           $(top_srcdir)/test/data/hess_bkg_lookup.fits

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

# Test data to be distributed in $(pyexecdir)/cscripts/tests/data/prod5-v0.1
testprod501data = $(top_srcdir)/test/data/prod5-v0.1/CTA-Performance-prod5-v0.1-South-20deg.FITS.tar

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
testprod501dir      = $(pyexecdir)/cscripts/tests/data/prod5-v0.1
testprod501_DATA    = $(testprod501data)

# Files to be distributed in the tar file
EXTRA_DIST = $(binscripts) $(pythonscripts) $(parfiles) \
             $(testscripts) $(testdata) $(testiactdata) \
             $(testiactevents) $(testprod501data)

# Add install hook rule to create symbolic links to cscripts without .py in the bindir
install-data-hook:
	@binlist='$(binscripts)';\
	for p in $$binlist; do \
		chmod +x "$(DESTDIR)$(bindir)/$$p"; \
	done

	@linklist1='$(linkscripts)';\
	for p in $$linklist1; do \
		rm -rf "$(DESTDIR)$(bindir)/$$p"; \
		ln -s "$(DESTDIR)$(pyexecdir)/cscripts/$$p.py" "$(DESTDIR)$(bindir)/$$p"; \
		chmod +x "$(DESTDIR)$(bindir)/$$p"; \
	done

# Uninstall hook rule to remove the cscripts from the bindir
uninstall-hook:
	@linklist2='$(linkscripts)';\
	for p in $$linklist2; do \
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
		rm -rf $(top_builddir)/cscripts/mputils.py; \
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

