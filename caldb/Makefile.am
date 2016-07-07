#############################################################################
# Makefile.am for calibration database                                      #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016 Juergen Knoedlseder                                    #
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

# Calibration files to install
caldb_prod2_indxdir      = $(datadir)/caldb/data/cta/prod2
caldb_prod2_indx_DATA    = $(srcdir)/data/cta/prod2/caldb.indx
#
caldb_prod2_north1dir    = $(datadir)/caldb/data/cta/prod2/bcf/North_0.5h
caldb_prod2_north1_DATA  = $(srcdir)/data/cta/prod2/bcf/North_0.5h/irf_file.fits.gz
caldb_prod2_north2dir    = $(datadir)/caldb/data/cta/prod2/bcf/North_5h
caldb_prod2_north2_DATA  = $(srcdir)/data/cta/prod2/bcf/North_5h/irf_file.fits.gz
caldb_prod2_north3dir    = $(datadir)/caldb/data/cta/prod2/bcf/North_50h
caldb_prod2_north3_DATA  = $(srcdir)/data/cta/prod2/bcf/North_50h/irf_file.fits.gz
#
caldb_prod2_south1dir    = $(datadir)/caldb/data/cta/prod2/bcf/South_0.5h
caldb_prod2_south1_DATA  = $(srcdir)/data/cta/prod2/bcf/South_0.5h/irf_file.fits.gz
caldb_prod2_south2dir    = $(datadir)/caldb/data/cta/prod2/bcf/South_5h
caldb_prod2_south2_DATA  = $(srcdir)/data/cta/prod2/bcf/South_5h/irf_file.fits.gz
caldb_prod2_south3dir    = $(datadir)/caldb/data/cta/prod2/bcf/South_50h
caldb_prod2_south3_DATA  = $(srcdir)/data/cta/prod2/bcf/South_50h/irf_file.fits.gz

# Calibration files to be distributed in the tar file
EXTRA_DIST = $(caldb_prod2_indx_DATA) \
             $(caldb_prod2_north1_DATA) \
             $(caldb_prod2_north2_DATA) \
             $(caldb_prod2_north3_DATA) \
             $(caldb_prod2_south1_DATA) \
             $(caldb_prod2_south2_DATA) \
             $(caldb_prod2_south3_DATA)
