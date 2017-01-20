#!/bin/bash -f
# ==========================================================================
# ctools Mac OS X package creation
#
# Copyright (C) 2017 Sylvie Brau-Nogué
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
CFITSIO=cfitsio3410
NCURSES=ncurses-6.0
READLINE=readline-7.0
GAMMALIB=gammalib-$VERSION
CTOOLS=ctools-$VERSION


# ============== #
# Set parameters #
# ============== #
WRKDIR=$PWD/osxpackage
INSTALLDIR=$WRKDIR/install
SRCDIR=$WRKDIR/src
PKGDIR=$WRKDIR/pkg
PRODDIR=$WRKDIR/prod
PLISTFILE=$PKGDIR/$CTOOLS-components.plist
DISTFILE=$PRODDIR/$CTOOLS.dist
DMGFILE=$WRKDIR/$CTOOLS-macosx10.7.dmg


# ====================== #
# Clean package creation #
# ====================== #
rm -rf $INSTALLDIR
rm -rf $PKGDIR
rm -rf $PRODDIR
rm -rf $DMGFILE


# ============================= #
# Create package directory tree #
# ============================= #
mkdir -p $SRCDIR
mkdir -p $PKGDIR
mkdir -p $PRODDIR


# Make sure that wget is installed
# TODO: If not installed, execute "brew install wget"


# =============== #
# Install ncurses #
# =============== #
cd $SRCDIR
if [ ! -d "$NCURSES" ]; then
    wget https://ftp.gnu.org/pub/gnu/ncurses/$NCURSES.tar.gz
    tar xvfz $NCURSES.tar.gz
fi
cd $NCURSES
./configure --prefix=$INSTALLDIR --with-shared
make
make install


# ================ #
# Install readline #
# ================ #
cd $SRCDIR
if [ ! -d "$READLINE" ]; then
  wget https://ftp.gnu.org/pub/gnu/readline/$READLINE.tar.gz
  tar xvfz $READLINE.tar.gz
fi
cd $READLINE
./configure --prefix=$INSTALLDIR
make
make install


# =============== #
# Install cfitsio #
# =============== #
cd $SRCDIR
if [ ! -d "cfitsio" ]; then
    wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/$CFITSIO.tar.gz
    tar xvfz $CFITSIO.tar.gz
fi
cd cfitsio
./configure --prefix=$INSTALLDIR
make shared
make install


# ================ #
# Install GammaLib #
# ================ #
cd $SRCDIR
if [ ! -d "$GAMMALIB" ]; then
    wget http://cta.irap.omp.eu/ctools/releases/gammalib/$GAMMALIB.tar.gz
    tar xvfz $GAMMALIB.tar.gz
fi
cd $GAMMALIB
./configure --prefix=$INSTALLDIR
make
make install


# ============== #
# Install ctools #
# ============== #
cd $SRCDIR
if [ ! -d "$CTOOLS" ]; then
    wget http://cta.irap.omp.eu/ctools/releases/ctools/$CTOOLS.tar.gz
    tar xvfz $CTOOLS.tar.gz
fi
cd $CTOOLS
./configure --prefix=$INSTALLDIR
make
make install


# ====================== #
# Build Mac OS X package #
# ====================== #
pkgbuild --analyze --root $INSTALLDIR $PLISTFILE
pkgbuild --identifier cta.irap.omp.eu.ctools.pkg \
         --version $VERSION \
         --root $INSTALLDIR \
         --install-location /usr/local/gamma \
         --component-plist $PLISTFILE \
         $PKGDIR/$CTOOLS.pkg


# ====================== #
# Build Mac OS X product #
# ====================== #
# Analyse existing package
productbuild --synthesize \
             --package $PKGDIR/$CTOOLS.pkg \
             $DISTFILE

# Add information to dist file
sed -i '' '$ i \
<license file="COPYING"/>\
<title>'$CTOOLS'</title>\
<organization>irap.omp.eu</organization>\
<description>CTA Science Analysis Tools</description>\
' "$PRODDIR/$CTOOLS.dist"

# Build product
productbuild --distribution $DISTFILE \
             --version $VERSION \
             --resources $SRCDIR/$CTOOLS \
             --package-path $WRKDIR \
             $PRODDIR/$CTOOLS.pkg

# Add additional files to production folder
cp $SRCDIR/$CTOOLS/COPYING $PRODDIR/License.txt

# Add ReadMe
/bin/cat <<EOM >$PRODDIR/ReadMe.txt
This package will install

                   ctools-$VERSION and gammalib-$VERSION

for Mac OS X 10.7 or later.

Installation requires approximately 260 MB of disk space. The package
will be installed in the directory /usr/local/gamma, hence you need
administrator privileges to do the install.

Simply click on the installer package and follow the instructions.
The installation will take less than one minute.

After installing the package, you have to initialise the ctools
environment. This will be done by invoking an initialisation script
that can be found in the /usr/local/gamma/bin directory.

To invoke this script you need to add the following to your .bashrc 
or $HOME/.profile script:

    export CTOOLS=/usr/local/gamma
    source \$CTOOLS/bin/gammalib-init.sh

If you use C shell or a variant then add the following to your
.cshrc or .tcshrc script:

    setenv CTOOLS /usr/local/gamma
    source \$CTOOLS/bin/gammalib-init.csh

More information on ctools can be found at
http://cta.irap.omp.eu/ctools/.

If you encounter problems in installing or using this package, please
send your error report to jurgen.knodlseder@irap.omp.eu. In your report,
please specify your Mac OS X version and Python version and the problem
you're facing.
EOM



# ==================== #
# Create Mac OS X disk #
# ==================== #
hdiutil create -volname $CTOOLS \
               -srcfolder $PRODDIR/$CTOOLS.pkg \
               -srcfolder $PRODDIR/License.txt \
               -srcfolder $PRODDIR/ReadMe.txt \
               $DMGFILE
