#!/bin/bash
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
# -------------------------------------------------------------------------
#
# This script installs the ctools and GammaLib packages in the location
# /usr/local/gamma and then builds a Mac OS X package which is finally put
# onto a disk image. Any existing content in /usr/local/gamma will be
# destroyed. The script needs the priviledges to create a folder in
# /usr/local.
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
NCURSES=ncurses-5.9
READLINE=readline-6.3
GAMMALIB=gammalib-$VERSION
CTOOLS=ctools-$VERSION


# ============== #
# Set parameters #
# ============== #
INSTALLDIR=/usr/local/gamma
WRKDIR=$PWD/pkg_build
SRCDIR=$WRKDIR/src
PKGDIR=$WRKDIR/pkg
PRODDIR=$WRKDIR/prod
PLISTFILE=$PKGDIR/$CTOOLS-components.plist
DISTFILE=$PRODDIR/$CTOOLS.dist
DMGFILE=$WRKDIR/$CTOOLS-macosx10.7.dmg
LOGFILE=$PWD/pkg_build.log


# ============================= #
# Secure installation directory #
# ============================= #
if [ -d "$INSTALLDIR" ]; then
    sudo mv $INSTALLDIR $INSTALLDIR.backup
fi


# ======================= #
# Clean package directory #
# ======================= #
sudo rm -rf $INSTALLDIR
rm -rf $PKGDIR
rm -rf $PRODDIR
rm -rf $DMGFILE
rm -rf $LOGFILE


# ============================= #
# Create package directory tree #
# ============================= #
mkdir -p $SRCDIR
mkdir -p $PKGDIR
mkdir -p $PRODDIR


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
        rm -rf $NCURSES.tar.gz
    fi
fi
cd $NCURSES
./configure --prefix=$INSTALLDIR --with-shared | tee -a $LOGFILE
make -j4 | tee -a $LOGFILE
sudo make install | tee -a $LOGFILE
make clean


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
        rm -rf $READLINE.tar.gz
    fi
fi
cd $READLINE
./configure --prefix=$INSTALLDIR | tee -a $LOGFILE
make -j4 | tee -a $LOGFILE
sudo make install | tee -a $LOGFILE
make clean


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
        rm -rf $CFITSIO.tar.gz
    fi
fi
cd cfitsio
./configure --prefix=$INSTALLDIR | tee -a $LOGFILE
make -j4 shared | tee -a $LOGFILE
sudo make install | tee -a $LOGFILE
make clean


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
        rm -rf $GAMMALIB.tar.gz
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
./configure --prefix=$INSTALLDIR | tee -a $LOGFILE
make -j4 | tee -a $LOGFILE
sudo make install | tee -a $LOGFILE
make clean


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
        rm -rf $CTOOLS.tar.gz
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
./configure --prefix=$INSTALLDIR | tee -a $LOGFILE
make -j4 | tee -a $LOGFILE
sudo make install | tee -a $LOGFILE
make clean


# ================= #
# Post-process code #
# ================= #
# Change links in libraries
for file in ${INSTALLDIR}/lib/*.dylib
do
    echo $file | tee -a $LOGFILE
    sudo install_name_tool -change "@rpath/libcfitsio.5.dylib" "$INSTALLDIR/lib/libcfitsio.5.dylib" $file
done

# Change links in binaries
for file in ${INSTALLDIR}/bin/ct*
do
    echo $file | tee -a $LOGFILE
    sudo install_name_tool -change "@rpath/libcfitsio.5.dylib" "$INSTALLDIR/lib/libcfitsio.5.dylib" $file
done

# Change links in Python modules
pythons="python2.7"
for python in $pythons
do
    # Print Python version
    echo $python | tee -a $LOGFILE
  
    # Change links in gammalib modules
    for file in ${INSTALLDIR}/lib/$python/site-packages/gammalib/_*.so
    do
        echo $file | tee -a $LOGFILE
        sudo install_name_tool -change "@rpath/libcfitsio.5.dylib" "$INSTALLDIR/lib/libcfitsio.5.dylib"  $file
    done

    # Change links in ctools modules
    for file in ${INSTALLDIR}/lib/$python/site-packages/ctools/_*.so
    do
        echo $file | tee -a $LOGFILE
        sudo install_name_tool -change "@rpath/libcfitsio.5.dylib" "$INSTALLDIR/lib/libcfitsio.5.dylib"  $file
    done
done


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
             --sign 'ctools' \
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


# ======================= #
# Clean package directory #
# ======================= #
rm -rf $PKGDIR
rm -rf $PRODDIR
