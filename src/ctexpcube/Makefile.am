#############################################################################
# Makefile.am for ctexpcube tool                                            #
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

# Set include directories
INCLUDES = -I$(top_srcdir)/src/support

# Program name
bin_PROGRAMS = ctexpcube

# Sources
sources = ctexpcube.cpp ctexpcube.hpp

# Libtool library
noinst_LTLIBRARIES = libctexpcube.la
libctexpcube_la_SOURCES = $(sources)
libctexpcube_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctexpcube_SOURCES = main.cpp
ctexpcube_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctexpcube_CPPFLAGS = @CPPFLAGS@
ctexpcube_LDADD = libctexpcube.la $(top_builddir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = $(srcdir)/ctexpcube.par

# Extra files to be distributed in the tar file
EXTRA_DIST = $(pfiles_DATA)
