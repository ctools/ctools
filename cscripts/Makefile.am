## Process this file with automake to produce Makefile.in

# Scripts to be distributed in $(pyexecdir)/cscripts
pythonscripts = __init__.py \
                cspull.py cssens.py cstsdist.py cscaldb.py csresmap.py \
                csiactobs.py csiactdata.py csobsdef.py csspec.py cslightcrv.py csobsinfo.py \
                csinfo.py csmodelinfo.py csmodelmerge.py cstsmapmerge.py \
                csobs2caldb.py \
                obsutils.py

# Scripts for which symbolic links will be created in bindir
linkscripts = cspull cssens cstsdist cscaldb csresmap \
              csiactobs csiactdata csobsdef csspec cslightcrv csobsinfo \
              csinfo csmodelinfo csmodelmerge cstsmapmerge \
              csobs2caldb

# Distribute scripts
EXTRA_DIST         = $(binscripts) $(pythonscripts)
binscriptsdir      = $(bindir)
binscripts_DATA    = $(binscripts)
pythonscriptsdir   = $(pyexecdir)/cscripts
pythonscripts_DATA = $(pythonscripts)

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
    
