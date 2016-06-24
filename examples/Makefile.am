..## Process this file with automake to produce Makefile.in

# Scripts to be distributed
scripts = README.md \
          make_pointings.py \
          make_spectrum.py \
          make_ts_distributions.py \
          pipeline_binned_disk.py \
          pipeline_binned_mem.py \
          pipeline_stacked_disk.py \
          pipeline_stacked_mem.py \
          pipeline_unbinned_disk.py \
          pipeline_unbinned_mem.py \
          show_butterfly.py \
          show_model.py \
          show_pha.py \
          show_pull_evolution.py \
          show_pull_histogram.py \
          show_response.py \
          show_sensitivity.py \
          show_spectrum.py \
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
    
