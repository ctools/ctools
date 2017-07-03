#############################################################################
# Makefile.am for ctphase tool                                              #
# ------------------------------------------------------------------------- #
# Copyright (C) 2017 Joshua Cardenzana                                      #
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

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctphase

# Sources
sources = ctphase.cpp ctphase.hpp

# Headers
pkginclude_HEADERS = ctphase.hpp

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctphase.la
libctphase_la_SOURCES = $(sources)
libctphase_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctphase_SOURCES = main.cpp
ctphase_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctphase_CPPFLAGS = @CPPFLAGS@
ctphase_LDADD = libctphase.la $(top_builddir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = $(srcdir)/ctphase.par

# Extra files to be distributed in the tar file
EXTRA_DIST = $(pfiles_DATA)
