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
WRKDIR=$PWD/pkgbuild
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
    if [ "$?" -ne "0" ]; then
        echo "*** Unable to download $NCURSES.tar.gz"
        exit 1
    else
        tar xvfz $NCURSES.tar.gz
    fi
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
    if [ "$?" -ne "0" ]; then
        echo "*** Unable to download $READLINE.tar.gz"
        exit 1
    else
        tar xvfz $READLINE.tar.gz
    fi
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
    if [ "$?" -ne "0" ]; then
        echo "*** Unable to download $CFITSIO.tar.gz"
        exit 1
    else
        tar xvfz $CFITSIO.tar.gz
    fi
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
    if [ "$?" -ne "0" ]; then
        if [ ! -d "gammalib" ]; then
            echo "*** Unable to download $GAMMALIB.tar.gz, try cloning from GitLab"
            export GIT_SSL_NO_VERIFY=true
            git clone https://cta-gitlab.irap.omp.eu/gammalib/gammalib.git
            if [ "$?" -ne "0" ]; then
                echo "*** Unable to clone GammaLib from GitLab"
                exit 1
            fi
        else
            echo "*** Unable to download $GAMMALIB.tar.gz, use code in gammalib directory"
        fi
        cd gammalib
    else
        tar xvfz $GAMMALIB.tar.gz
        cd $GAMMALIB
    fi
fi
# If we got code from GitLab, then try using release branch, and if not
# available, use devel branch.
USE_BRANCH=
if [ -d ".git" ]; then
    git fetch
    git checkout release
    if [ "$?" -ne "0" ]; then
        echo "*** No release branch found on GitLab, try devel branch"
        git checkout devel
        if [ "$?" -ne "0" ]; then
            echo "*** No devel branch found on GitLab"
            exit 1
        else
            echo "*** Use devel branch of GammaLib"
            USE_BRANCH=devel
        fi
    else
        echo "*** Use release branch of GammaLib"
        USE_BRANCH=release
    fi
    git pull
    ./autogen.sh
fi
# Build code
./configure --prefix=$INSTALLDIR
make
make install


# ============== #
# Install ctools #
# ============== #
cd $SRCDIR
CTOOLS_DIR=$CTOOLS
if [ ! -d "$CTOOLS" ]; then
    wget http://cta.irap.omp.eu/ctools/releases/ctools/$CTOOLS.tar.gz
    if [ "$?" -ne "0" ]; then
        if [ ! -d "ctools" ]; then
            echo "*** Unable to download $CTOOLS.tar.gz, try cloning from GitLab"
            export GIT_SSL_NO_VERIFY=true
            git clone https://cta-gitlab.irap.omp.eu/ctools/ctools.git
            if [ "$?" -ne "0" ]; then
                echo "*** Unable to clone ctools from GitLab"
                exit 1
            fi
        else
            echo "*** Unable to download $CTOOLS.tar.gz, use code in ctools directory"
        fi
        cd ctools
        CTOOLS_DIR=ctools
    else
        tar xvfz $CTOOLS.tar.gz
        cd $CTOOLS
    fi
fi
# If we got code from GitLab, then use the same branch that was used for
# GammaLib
if [ -d ".git" ]; then
    git fetch
    git checkout $USE_BRANCH
    if [ "$?" -ne "0" ]; then
        echo "*** No $USE_BRANCH branch found on GitLab"
        exit 1
    fi
    git pull
    ./autogen.sh
fi
# Build code
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
             --package-path $PKGDIR \
             $PRODDIR/$CTOOLS.pkg

# Add additional files to production folder
cp $SRCDIR/$CTOOLS_DIR/COPYING $PRODDIR/License.txt

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
