## Process this file with automake to produce Makefile.in

# Scripts to be distributed in bindir
binscripts = cssens

# Scripts to be distributed in pythondir
pythonscripts = obsutils.py

# Distribute scripts
EXTRA_DIST         = $(binscripts) $(pythonscripts)
binscriptsdir      = $(bindir)
binscripts_DATA    = $(binscripts)
pythonscriptsdir   = $(pythondir)
pythonscripts_DATA = $(pythonscripts)

# Add hook rule to make binscripts executable
install-data-hook:
	@list='$(binscripts_DATA)';\
	for p in $$list; do \
	  if test -f "$$p"; then chmod +x "$(DESTDIR)$(binscriptsdir)/$$p"; fi; \
	done
