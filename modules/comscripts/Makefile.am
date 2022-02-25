#############################################################################
# Makefile for comscripts Python scripts                                    #
# ------------------------------------------------------------------------- #
# Copyright (C) 2021-2022 Juergen Knoedlseder                               #
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

# Scripts to be distributed in $(pyexecdir)/comscripts
pythonscripts = $(top_builddir)/modules/comscripts/__init__.py \
                $(srcdir)/comgendb.py \
                $(srcdir)/comlixfit.py \
                $(srcdir)/comlixmap.py \
                $(srcdir)/comobsadd.py \
                $(srcdir)/comobsback.py \
                $(srcdir)/comobsbin.py \
                $(srcdir)/comobsmodel.py \
                $(srcdir)/comobsres.py \
                $(srcdir)/comobsselect.py \
                $(srcdir)/comobssim.py \
                $(srcdir)/compulbin.py \
                $(srcdir)/comsrcdetect.py

# Parfiles to be distributed in $(prefix)/syspfiles
parfiles = $(srcdir)/comgendb.par \
           $(srcdir)/comlixfit.par \
           $(srcdir)/comlixmap.par \
           $(srcdir)/comobsadd.par \
           $(srcdir)/comobsback.par \
           $(srcdir)/comobsbin.par \
           $(srcdir)/comobsmodel.par \
           $(srcdir)/comobsres.par \
           $(srcdir)/comobsselect.par \
           $(srcdir)/comobssim.par \
           $(srcdir)/compulbin.par \
           $(srcdir)/comsrcdetect.par

# Scripts for which symbolic links will be created in bindir
linkscripts = comgendb \
              comlixfit \
              comlixmap \
              comobsadd \
              comobsback \
              comobsbin \
              comobsmodel \
              comobsres \
              comobsselect \
              comobssim \
              compulbin \
              comsrcdetect

# Test scripts to be distributed in $(pyexecdir)/comscripts/tests
testscripts = $(srcdir)/tests/__init__.py \
              $(top_srcdir)/test/testing.py \
              $(top_srcdir)/test/test_python_comscripts.py \
              $(top_srcdir)/test/test_comlixfit.py \
              $(top_srcdir)/test/test_comlixmap.py \
              $(top_srcdir)/test/test_comobsadd.py \
              $(top_srcdir)/test/test_comobsback.py \
              $(top_srcdir)/test/test_comobsbin.py \
              $(top_srcdir)/test/test_comobsmodel.py \
              $(top_srcdir)/test/test_comobsres.py \
              $(top_srcdir)/test/test_comobsselect.py \
              $(top_srcdir)/test/test_comobssim.py \
              $(top_srcdir)/test/test_comsrcdetect.py

# Test data to be distributed in $(pyexecdir)/comscripts/tests/data/comptel
testdata = $(top_srcdir)/test/data/comptel/ebounds_std.fits.gz \
           $(top_srcdir)/test/data/comptel/ebounds_std1.fits.gz \
           $(top_srcdir)/test/data/comptel/models.xml \
           $(top_srcdir)/test/data/comptel/obs_binned.xml \
           $(top_srcdir)/test/data/comptel/obs_binned1.xml \
           $(top_srcdir)/test/data/comptel/obs_selected.xml \
           $(top_srcdir)/test/data/comptel/tsmap.fits.gz

# Test datastore to be distributed in $(pyexecdir)/comscripts/tests/data/comptel/datastore
testdatastore = $(top_srcdir)/test/data/comptel/datastore/iaq_000750-001000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/iaq_001000-003000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/iaq_003000-010000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/iaq_010000-030000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_drb_psd0-110_000750-001000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_drb_psd0-110_001000-003000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_drb_psd0-110_003000-010000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_drb_psd0-110_010000-030000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_dre_psd0-110_000750-001000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_dre_psd0-110_001000-003000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_dre_psd0-110_003000-010000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_dre_psd0-110_010000-030000keV.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_drg.fits.gz \
                $(top_srcdir)/test/data/comptel/datastore/vp0001_0_drx.fits.gz

# Test data to be distributed in $(pyexecdir)/comscripts/tests/data/comptel/data/dbase
testdatabase = $(top_srcdir)/test/data/comptel/data/dbase/dbase.fits.gz

# Test data to be distributed in $(pyexecdir)/comscripts/tests/data/comptel/data/dbase/tim
testdatatim = $(top_srcdir)/test/data/comptel/data/dbase/tim/vp0001_0_tim.fits.gz

# Test data to be distributed in $(pyexecdir)/comscripts/tests/data/comptel/data/dbase/xml
testdataxml = $(top_srcdir)/test/data/comptel/data/dbase/xml/vp0001_0.xml

# Test data to be distributed in $(pyexecdir)/comscripts/tests/data/comptel/data/phase01/vp0001_0
testdatavp = $(top_srcdir)/test/data/comptel/data/phase01/vp0001_0/vp0001_evp.fits.gz \
             $(top_srcdir)/test/data/comptel/data/phase01/vp0001_0/vp0001_oad.fits.gz

# Install scripts, parfiles, test scripts and test data. Note that we don't
# need to prepend $(DESTDIR) here (not fully clear why, but that's the way
# it is)
binscriptsdir       = $(bindir)
binscripts_DATA     = $(binscripts)
pythonscriptsdir    = $(pyexecdir)/comscripts
pythonscripts_DATA  = $(pythonscripts)
pfilesdir           = $(prefix)/syspfiles
pfiles_DATA         = $(parfiles)
testscriptsdir      = $(pyexecdir)/comscripts/tests
testscripts_DATA    = $(testscripts)
testdatadir         = $(pyexecdir)/comscripts/tests/data/comptel
testdata_DATA       = $(testdata)
testdatastoredir    = $(pyexecdir)/comscripts/tests/data/comptel/datastore
testdatastore_DATA  = $(testdatastore)
testdatabasedir     = $(pyexecdir)/comscripts/tests/data/comptel/data/dbase
testdatabase_DATA   = $(testdatabase)
testdatatimdir      = $(pyexecdir)/comscripts/tests/data/comptel/data/dbase/tim
testdatatim_DATA    = $(testdatatim)
testdataxmldir      = $(pyexecdir)/comscripts/tests/data/comptel/data/dbase/xml
testdataxml_DATA    = $(testdataxml)
testdatavpdir       = $(pyexecdir)/comscripts/tests/data/comptel/data/phase01/vp0001_0
testdatavp_DATA     = $(testdatavp)

# Files to be distributed in the tar file
EXTRA_DIST = $(binscripts) $(pythonscripts) $(parfiles) \
             $(testscripts) $(testdata) $(testdatastore) $(testdatabase) \
             $(testdatatim) $(testdataxml) $(testdatavp)

# Add install hook rule to create symbolic links to comscripts without .py in the bindir
install-data-hook:
	@binlist='$(binscripts)';\
	for p in $$binlist; do \
		chmod +x "$(DESTDIR)$(bindir)/$$p"; \
	done

	@linklist1='$(linkscripts)';\
	for p in $$linklist1; do \
		rm -rf "$(DESTDIR)$(bindir)/$$p"; \
		ln -s "$(DESTDIR)$(pyexecdir)/comscripts/$$p.py" "$(DESTDIR)$(bindir)/$$p"; \
		chmod +x "$(DESTDIR)$(bindir)/$$p"; \
	done

# Uninstall hook rule to remove the comscripts from the bindir
uninstall-hook:
	@linklist2='$(linkscripts)';\
	for p in $$linklist2; do \
		rm -rf "$(DESTDIR)$(bindir)/$$p"; \
	done

# Files to be removed in build directory
clean-local:
	@if ! [ '@top_srcdir@' -ef '@top_builddir@' ]; then \
		rm -rf $(top_builddir)/modules/comscripts/com*.py; \
		rm -rf $(top_builddir)/modules/comscripts/*.par; \
	fi
	rm -rf $(top_builddir)/modules/comscripts/*.pyc
	rm -rf $(top_builddir)/modules/comscripts/__pycache__

# Uninstall the comscripts Python module
uninstall-local:
	rm -rf $(DESTDIR)$(pyexecdir)/comscripts/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/comscripts/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/comscripts/__pycache__
	rm -rf $(DESTDIR)$(pyexecdir)/comscripts/tests/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/comscripts/tests/__pycache__

