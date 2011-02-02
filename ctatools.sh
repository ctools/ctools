#! /bin/bash -f

sourcedir=ctatools
installdir=$(pwd)/build

echo "Remove existing folder"
rm -rf $sourcedir
#rm -rf $installdir

echo "Checkout latest ctatools"
cvs -Q -d /home/cvs export -N -r HEAD "ctatools"

echo "Create Makefile.in and configure scripts"
cd $sourcedir
#
mkdir m4
aclocal -I m4
libtoolize --copy
autoconf
autoheader
automake --add-missing --copy
#
rm -f autogen.sh      
#rm -rf m4
rm -rf autom4te.cache
#rm -f configure.ac
#rm -f Makefile.am
#
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
export PKG_CONFIG_PATH=/home/knodlseder/export/gammalib_build/lib/pkgconfig
./configure --prefix=$installdir

echo "Compile ctatools"
make -j4

echo "Install ctatools"
make install

#echo "Check GammaLib"
#export PATH=$installdir/bin:$PATH
#export LD_LIBRARY_PATH=$installdir/lib:$LD_LIBRARY_PATH
#export PYTHONPATH=$installdir/lib/python2.5/site-packages:$PYTHONPATH
#make check
