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
    mv $INSTALLDIR $INSTALLDIR.backup
fi


# ====================== #
# Clean package creation #
# ====================== #
#umount $WRKDIR
#rm -rf $WRKDIR


# ============================= #
# Create package directory tree #
# ============================= #
#mkdir -p $WRKDIR


# ================================= #
# Create Mac OS X RAM disk (488 MB) #
# ================================= #
#DEVICE=$(hdiutil attach ram://1000000 -nomount)
#diskutil erasevolume HFS+ 'ctools-test' $DEVICE
#umount -f /Volumes/ctools-test
#sudo diskutil enableOwnership $DEVICE
#mount -t hfs $DEVICE $WRKDIR


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
python -c 'import gammalib; gammalib.test()'
python -c 'import ctools; ctools.test()'
python -c 'import cscripts; cscripts.test()'


# ======================= #
# Clean package directory #
# ======================= #
#rm -rf $INSTALLDIR


# ============================== #
# Recover installation directory #
# ============================== #
#if [ -d "$INSTALLDIR.backup" ]; then
#    mv $INSTALLDIR.backup $INSTALLDIR
#fi


# =============== #
# Detach RAM disk #
# =============== #
#hdiutil detach $WRKDIR

