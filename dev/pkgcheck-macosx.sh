#!/bin/bash
# ==========================================================================
# ctools Mac OS X package validation
#
# Copyright (C) 2017 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------
#
# This script checks the ctools Mac OS X package by installing it into the
# /usr/local/gamma directory and running the Python tests.
#
# ==========================================================================

# ============================================= #
# Get ctools/GammaLib version from command line #
# ============================================= #
if [ $# -ne 1 ] ;  then
    echo "Please specify ctools/GammaLib version in the form x.y.z"
    exit 1
fi
VERSION=$1


# =============================== #
# Set software component versions #
# =============================== #
CTOOLS=ctools-$VERSION


# ============== #
# Set parameters #
# ============== #
INSTALLDIR=/usr/local/gamma
WRKDIR=$PWD/pkg_check
PKGDIR=$PWD/pkg_build
DMGFILE=$PKGDIR/$CTOOLS-macosx10.7.dmg
LOGFILE=$PWD/pkg_check.log


# ============================= #
# Secure installation directory #
# ============================= #
if [ -d "$INSTALLDIR" ]; then
    sudo mv $INSTALLDIR $INSTALLDIR.backup
fi


# ======================= #
# Clean working directory #
# ======================= #
rm -rf $WRKDIR


# ====================================== #
# Create and step into working directory #
# ====================================== #
mkdir -p $WRKDIR
cd $WRKDIR


# =============== #
# Install package #
# =============== #
hdiutil attach $DMGFILE
sudo installer -pkg /Volumes/$CTOOLS/$CTOOLS.pkg -target /
hdiutil detach /Volumes/$CTOOLS


# ================= #
# Configure package #
# ================= #
export GAMMALIB=$INSTALLDIR
source $GAMMALIB/bin/gammalib-init.sh
export CTOOLS=$INSTALLDIR
source $CTOOLS/bin/ctools-init.sh


# ============ #
# Test package #
# ============ #
python -c 'import gammalib; gammalib.test()' | tee -a $LOGFILE
python -c 'import ctools; ctools.test()' | tee -a $LOGFILE
python -c 'import cscripts; cscripts.test()' | tee -a $LOGFILE


# ======================= #
# Clean package directory #
# ======================= #
sudo rm -rf $INSTALLDIR


# ============================== #
# Recover installation directory #
# ============================== #
if [ -d "$INSTALLDIR.backup" ]; then
    sudo mv $INSTALLDIR.backup $INSTALLDIR
fi
