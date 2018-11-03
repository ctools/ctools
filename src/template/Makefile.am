#############################################################################
# Makefile.am for template module                                           #
# ------------------------------------------------------------------------- #
# Copyright (C) 2018 Juergen Knoedlseder                                    #
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
bin_PROGRAMS = ctool_base \
               ctool_obs \
               ctool_like

# Program definition
ctool_base_SOURCES = main_base.cpp ctool_base.cpp
ctool_base_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@ @OPENMP_LIB_LDFLAGS@
ctool_base_CPPFLAGS = @CPPFLAGS@
ctool_base_LDADD = $(top_builddir)/src/support/libsupport.la

# Program definition
ctool_obs_SOURCES = main_obs.cpp ctool_obs.cpp
ctool_obs_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@ @OPENMP_LIB_LDFLAGS@
ctool_obs_CPPFLAGS = @CPPFLAGS@
ctool_obs_LDADD = $(top_builddir)/src/support/libsupport.la

# Program definition
ctool_like_SOURCES = main_like.cpp ctool_like.cpp
ctool_like_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@ @OPENMP_LIB_LDFLAGS@
ctool_like_CPPFLAGS = @CPPFLAGS@
ctool_like_LDADD = $(top_builddir)/src/support/libsupport.la

# Extra files to be distributed in the tar file
EXTRA_DIST = cscript.par \
             cscript.rst \
             cscript_base.py \
             cscript_like.py \
             cscript_obs.py \
             cscript_test.py \
             ctool.par \
             ctool.rst \
             ctool_base.cpp \
             ctool_base.hpp \
             ctool_base.i \
             ctool_like.cpp \
             ctool_like.hpp \
             ctool_like.i \
             ctool_Makefile.am \
             ctool_obs.cpp \
             ctool_obs.hpp \
             ctool_obs.i \
             ctool_test.py \
             main_base.cpp \
             main_like.cpp \
             main_obs.cpp
