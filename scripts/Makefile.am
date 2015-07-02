## Process this file with automake to produce Makefile.in

# Scripts to be distributed in bindir
binscripts = cta_kb_bkgrates.py cta_root2irf.py cta_make_nodes.py \
             tsmerge.py

# Scripts to be distributed in pythondir
pythonscripts = cspull.py cssens.py cstsdist.py cscaldb.py csresmap.py cshessobs.py \
                csobsdef.py csspec.py cslightcrv.py \
                lsmodel.py

# Scripts for which symbolic links will be created in bindir
linkscripts = cspull cssens cstsdist cscaldb csresmap \
              csobsdef csspec cshessobs cslightcrv \
              lsmodel csinfo

# Distribute scripts
EXTRA_DIST         = $(binscripts) $(pythonscripts)
binscriptsdir      = $(bindir)
binscripts_DATA    = $(binscripts)
pythonscriptsdir   = $(pythondir)/ctools
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
      ln -s "$(DESTDIR)$(pythondir)/ctools/$$p.py" "$(DESTDIR)$(bindir)/$$p"; \
	  chmod +x "$(DESTDIR)$(pythondir)/ctools/$$p.py"; \
	done
    
