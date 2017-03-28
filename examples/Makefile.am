#############################################################################
# Makefile for ctools examples                                              #
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

# Scripts to be install
scripts = README.md \
          generate_prod3_irfs.py \
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
          show_irf.py \
          show_model.py \
          show_pha.py \
          show_pull_evolution.py \
          show_pull_histogram.py \
          show_response.py \
          show_sensitivity.py \
          show_spectrum.py \
          show_lightcurve.py \
          show_ts_distribution.py \
          show_obs.py \
          show_pointings.py

# Set install directories
scriptsdir    = $(datadir)/examples/python
scripts_DATA  = $(scripts)

# Files to be distributed in the tar file
EXTRA_DIST = $(scripts)

# Add hook rule to make the examples executable
install-data-hook:
	@list='$(scripts)';\
	for p in $$list; do \
	  chmod +x "$(DESTDIR)$(scriptsdir)/$$p"; \
	done
    
