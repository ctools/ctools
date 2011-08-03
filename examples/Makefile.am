## Process this file with automake to produce Makefile.in

# Scripts to be distributed
scripts = README \
          make_binned_analysis.py \
          make_pull_at_sensitivity_limit.py \
          make_ts_distributions.py \
          make_unbinned_analysis.py \
          show_pull_histogram.py \
          show_ts_distribution.py

# Set directories
scriptsdir   = $(datadir)/examples/python
scripts_DATA = $(scripts)

# Extra files to be distributed in the tar file
EXTRA_DIST = $(scripts)

# Add hook rule to make the examples executable
install-data-hook:
	@list='$(scripts)';\
	for p in $$list; do \
	  chmod +x "$(DESTDIR)$(scriptsdir)/$$p"; \
	done
    
