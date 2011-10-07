#! /bin/bash -f

# Set parameters
installdir=$(pwd)/build

echo "Checkout latest ctools"
rm -rf ctools
cvs -Q -d /home/cvs export -N -r HEAD "ctools"

# Extract version number
version=`cat ctools/configure.ac | grep 'AC_INIT' | awk -F"[" '{print $3}' | sed 's/],//' | sed 's/\./ /g'`
version=`printf "%2.2d-%2.2d-%2.2d" $version`
echo "ctools version: "$version

# Set directory name
sourcedir=ctools-$version
rm -rf $sourcedir
mv ctools $sourcedir

echo "Create Makefile.in and configure scripts"
cd $sourcedir
./autogen.sh
rm -f ctools.sh
cd ..

echo "Create ctool_wrap.cpp and ctool.py files for python binding (making swig obsolete)"
swig -c++ -python -Wall -includeall -I$installdir/share/gammalib/swig -o $sourcedir/pyext/ctools_wrap.cpp -outdir $sourcedir/pyext $sourcedir/pyext/ctools.i

echo "Make files read/write"
chmod -R u+rw $sourcedir

echo "Create tarball"
tar cvf - $sourcedir > $sourcedir.tar
rm -f $sourcedir.tar.gz
gzip $sourcedir.tar

echo "Configure ctools"
cd $sourcedir
./configure --prefix=$installdir

echo "Compile ctools"
make -j10

echo "Check ctools"
make check

echo "Install ctools"
make install
