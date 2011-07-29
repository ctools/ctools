## Process this file with automake to produce Makefile.in

# Scripts for which symbolic links should be created in bindir
binscripts = cspull cssens

# Scripts to be distributed in pythondir
pythonscripts = cspull.py cssens.py obsutils.py

# Distribute scripts
EXTRA_DIST         = $(pythonscripts)
pythonscriptsdir   = $(pythondir)
pythonscripts_DATA = $(pythonscripts)

# Add hook rule to create symbolic links for cscripts
install-data-hook:
	@list='$(binscripts)';\
	for p in $$list; do \
      ln -s "$(DESTDIR)$(pythondir)/$$p.py" "$(DESTDIR)$(bindir)/$$p"; \
	  if test -f "$$p"; then chmod +x "$(DESTDIR)$(binscriptsdir)/$$p"; fi; \
	  chmod +x "$(DESTDIR)$(pythondir)/$$p.py"; \
	done
