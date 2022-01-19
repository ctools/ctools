#############################################################################
# Makefile for models                                                       #
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

# Models to install
models = 1dc_howto.xml \
         1dc_howto_unbinned.xml \
         1dc_howto_vela.xml \
         1dc_howto_velapsr.xml \
         comptel_howto_gc.xml \
         bkg_irf.xml \
         bkg_cube.xml \
         test_irf.xml \
         test_cube.xml \
         crab.xml \
         crab_beyond.xml \
         crab_2bkg.xml \
         crab_nobkg.xml \
         crab_poly.xml \
         crab_profile.xml \
         disk.xml \
         gauss.xml \
         shell.xml \
         hess1825_26.xml \
         hess_bkg_lookup.fits \
         ephem_xte.psrtime \
         README.md

# Install models
modelsdir   = $(datadir)/models
models_DATA = $(models)

# Extra files to be distributed in the tar file
EXTRA_DIST = $(models)
