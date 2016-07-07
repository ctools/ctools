#############################################################################
# Makefile.am for ctulimit tool                                             #
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
bin_PROGRAMS = ctulimit

# Sources
sources = ctulimit.cpp ctulimit.hpp

# Libtool library
noinst_LTLIBRARIES = libctulimit.la
libctulimit_la_SOURCES = $(sources)
libctulimit_la_LDFLAGS = @OPENMP_LIB_LDFLAGS@

# Program definition
ctulimit_SOURCES = main.cpp
ctulimit_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctulimit_CPPFLAGS = @CPPFLAGS@
ctulimit_LDADD = libctulimit.la $(top_builddir)/src/support/libsupport.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = $(srcdir)/ctulimit.par

# Extra files to be distributed in the tar file
EXTRA_DIST = $(pfiles_DATA)
