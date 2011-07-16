#! /bin/bash -f

# Set parameters
installdir=$(pwd)/build

echo "Checkout latest ctatools"
rm -rf ctatools
cvs -Q -d /home/cvs export -N -r HEAD "ctatools"

# Extract version number
version=`cat ctatools/configure.ac | grep 'AC_INIT' | awk -F"[" '{print $3}' | sed 's/],//' | sed 's/\./ /g'`
version=`printf "%2.2d-%2.2d-%2.2d" $version`
echo "ctatools version: "$version

# Set directory name
sourcedir=ctatools-$version
rm -rf $sourcedir
mv ctatools $sourcedir

echo "Create Makefile.in and configure scripts"
cd $sourcedir
mkdir m4
./autogen.sh
#rm -rf autom4te.cache
rm -f ctatools.sh
cd ..

echo "Create ctatool_wrap.cpp and ctatool.py files for python binding (making swig obsolete)"
swig -c++ -python -Wall -includeall -I$installdir/share/gammalib/swig -o $sourcedir/pyext/ctatools_wrap.cpp -outdir $sourcedir/pyext $sourcedir/pyext/ctatools.i

echo "Make files read/write"
chmod -R u+rw $sourcedir

echo "Create tarball"
tar cvf - $sourcedir > $sourcedir.tar
rm -f $sourcedir.tar.gz
gzip $sourcedir.tar

echo "Configure ctatools"
cd $sourcedir
./configure --prefix=$installdir

echo "Compile ctatools"
make -j10

echo "Check ctatools"
make check

echo "Install ctatools"
make install
