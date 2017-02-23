#!/bin/bash -f
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
#CFITSIO=cfitsio3410
#NCURSES=ncurses-6.0
#READLINE=readline-7.0
#GAMMALIB=gammalib-$VERSION
CTOOLS=ctools-$VERSION


# ============== #
# Set parameters #
# ============== #
WRKDIR=$PWD/pkg_check
PKGDIR=$PWD/pkg_build


#INSTALLDIR=$WRKDIR/install
#SRCDIR=$WRKDIR/src
#PKGDIR=$WRKDIR/pkg
#PRODDIR=$WRKDIR/prod
#PLISTFILE=$PKGDIR/$CTOOLS-components.plist
#DISTFILE=$PRODDIR/$CTOOLS.dist
DMGFILE=$PKGDIR/$CTOOLS-macosx10.7.dmg


# ====================== #
# Clean package creation #
# ====================== #
umount $WRKDIR
rm -rf $WRKDIR


# ============================= #
# Create package directory tree #
# ============================= #
mkdir -p $WRKDIR


# ================================= #
# Create Mac OS X RAM disk (488 MB) #
# ================================= #
DEVICE=$(hdiutil attach ram://1000000 -nomount)
diskutil erasevolume HFS+ 'ctools-test' $DEVICE
umount -f /Volumes/ctools-test
#sudo diskutil enableOwnership $DEVICE       # You must be root for this command
mount -t hfs $DEVICE $WRKDIR


# =============== #
# Install package #
# =============== #
echo $DMGFILE
hdiutil attach $DMGFILE
sudo installer -pkg /Volumes/$CTOOLS/$CTOOLS.pkg -target $WRKDIR # You must be root for this command
hdiutil detach /Volumes/$CTOOLS


# ================= #
# Configure package #
# ================= #
env
export GAMMALIB=$WRKDIR/usr/local/gamma
source $GAMMALIB/bin/gammalib-init.sh
export CTOOLS=$WRKDIR/usr/local/gamma
source $CTOOLS/bin/ctools-init.sh
echo
echo
env


# ============ #
# Test package #
# ============ #
python -c 'import gammalib; gammalib.test()'
python -c 'import ctools; ctools.test()'


# =============== #
# Detach RAM disk #
# =============== #
#hdiutil detach $WRKDIR

