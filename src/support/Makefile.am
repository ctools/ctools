#############################################################################
# Makefile.am for ctool base classes                                        #
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

# Define sources for this directory
sources = support.cpp support.hpp \
          ctool.cpp ctool.hpp \
          ctobservation.cpp ctobservation.hpp \
          ctlikelihood.cpp ctlikelihood.hpp

# Headers
pkginclude_HEADERS = support.hpp \
                     ctool.hpp \
                     ctobservation.hpp \
                     ctlikelihood.hpp

# Build libtool library
noinst_LTLIBRARIES = libsupport.la
libsupport_la_SOURCES = $(sources)
libsupport_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@
