ctools information
==================
* Version:             1.3.0.dev1 (12 August 2016)
* Author:              Juergen Knoedlseder (jurgen.knodlseder@irap.omp.eu)
* GammaLib dependency: 1.3.0.dev1

[![Build Status](https://cta-jenkins.irap.omp.eu/buildStatus/icon?job=ctools-integrate-os)](https://cta-jenkins.irap.omp.eu/job/ctools-integrate-os/)

[![Quality Gate](https://cta-sonar.irap.omp.eu/api/badges/gate?key=ctools)](https://cta-sonar.irap.omp.eu/dashboard/index/ctools)


License information
===================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


What's new in this release?
===========================
See the files [NEWS](NEWS) and [ChangeLog](ChangeLog).


What are the ctools anyway?
===========================
ctools are ftools-like executable for the scientific analysis of
CTA observations.  They are based on GammaLib, a versatile toolbox 
for the high-level analysis of astronomical gamma-ray data.

The following tools are available:

    ctbin       - event binning
    ctbkgcube   - generate a background cube
    ctbutterfly - create a butterfly
    ctcubemask  - mask bins in binned analysis
    ctedispcube - generate energy dispersion cube
    cterror     - likelihood profile error estimation
    ctexpcube   - generate an exposure cube
    ctlike      - maximum likelihood model fitting
    ctmapcube   - generate sky map cube
    ctmodel     - generation of model counts map
    ctobssim    - simulation of CTA observations
    ctpsfcube   - generate a PSF cube
    ctselect    - event selection
    ctskymap    - CTA sky mapping tool
    cttsmap     - generate a TS map
    ctulimit    - compute upper limits


Web sites
=========
http://cta.irap.omp.eu/ctools                     - for ctools users
https://cta-redmine.irap.omp.eu/projects/ctools   - for ctools development


Prerequisites
=============
ctools require GammaLib.  Please refer to http://gammalib.sourceforge.net
for instructions about how to install GammaLib.  The GammaLib code can
be downloaded from https://sourceforge.net/projects/gammalib/.

Once GammaLib is properly installed, make sure that you added the setup
script to your .bashrc or $HOME/.profile script:

    export GAMMALIB=/usr/local/gamma
    source $GAMMALIB/bin/gammalib-init.sh

If you use C shell or a variant then add the following to your
.cshrc or .tcshrc script:

    setenv GAMMALIB /usr/local/gamma
    source $GAMMALIB/bin/gammalib-init.csh

If you have installed GammaLib in another directory than /usr/local/gamma,
please adapt the path correspondingly.

If you really insist, you may install ctools in a directory different to
that hosting GammaLib, but we highly recommend to install both packages
together.


Unix Installation
=================
To build and install ctools, simply type the following:

     $ ./configure
     $ make
     $ make install

By default ctools installs itself in /usr/local/gamma.  If you need to
install ctools in a different location or in your home directory, use
the --prefix option to ./configure.  For example:

     $ ./configure --prefix=/home/yourname/projects
     $ make
     $ make install

The file INSTALL details more about using configure. Also try

     $ ./configure --help.


Macintosh OS X Installation
============================
ctools is known to work on various flavors of OS X.  To cope with
different system versions and architectures, there are two Mac
specific configure options:

     $ ./configure --enable-universalsdk[=PATH]

creates a universal build of ctools.  The optional argument specifies
which OSX SDK should be used to perform the build.  This defaults to
"/Developer/SDKs/MacOSX.10.4u.sdk".  Specify "/" when building on a 10.5
system or higher, especially when building 64-bit code.

     $ ./configure --with-univeral-archs=VALUE

specifies the kind of universal build that should be created.  Possible
values are: "32-bit", "3-way", "intel" or "all".  By default, a "32-bit" 
build will be made.  This option is only valid when
"--enable-universalsdk" is specified.

These options are in particular needed if your Python architecture differs
from the default architecture of your system.  To examine the Python
architecture you may type:

     $ file `which python`

which will return the architectures that are compiled in the Mach-0
executable:

     i386    32-bit intel
     ppc     32-bit powerpc
     ppc64   64-bit powerpc
     x86_64  64-bit intel

If Python is 32-bit (ppc, i386) but the compiler produces by default
64-bit code (ppc64, x86_64), the Python module will not work.  Using

     $ ./configure --enable-universalsdk=/

will force a universal 32-bit build which creates code for ppc and 
i386.  If on the other hand Python is 64-bit (ppc64, x86_64) but the
compiler produces by default 32-bit code (ppc, i386), the option

     $ ./configure --enable-universalsdk=/ --with-univeral-archs=3-way

will generate a universal build which contains 32-bit and 64-bit code.


BSD Installation
================
ctools has been tested on FreeBSD successfully.  Follow the Linux
installation instructions above.


Solaris Installation
====================
ctools compile on Solairs, but there is an issue with using it as a
shared library (see "Known problems" below).


Windows Installation
====================
There have been no efforts so far to compile ctools under Windows.


Testing
=======
If you want to test ctools before installation, type the following:

     $ make check


Setting up your environment
===========================
Before using ctools you have to setup some environment variables.
This will be done automatically by an initialisation script that will
be installed in the bin directory.

Assuming that you have installed ctools in the default directory 
`/usr/local/gamma` you need to add the following to your `$HOME/.bashrc` or 
`$HOME/.profile` script on a Linux machine:

    export CTOOLS=/usr/local/gamma
    source $CTOOLS/bin/ctools-init.sh

If you use C shell or a variant then add the following to your
`$HOME/.cshrc` or `$HOME/.tcshrc` script:

    setenv CTOOLS /usr/local/gamma
    source $CTOOLS/bin/ctools-init.csh


Getting started
===============
See the online documentation at http://cta.irap.omp.eu/ctools.


Documentation
=============
The doc directory (usually at /usr/local/gamma/share/doc/ctools)  
contains the most recent set of updated documentation for this release.
A detailed documentation can be created by typing:

     $ make doc

before installing the library. Two types of documentation exist:
* code documentation
* user documentation

Code documentation is created using Doxygen.  You need Doxygen on your
system to generate the code documentation.  This includes man pages.  
Doxygen can be obtained from http://www.stack.nl/~dimitri/doxygen/. 

User documentation is created using Sphinx.  You need Sphinx on your
system to generate the user documentation.  Sphinx can be obtained
from http://sphinx-doc.org/install.html.


Bug reports
===========
To report or search for bugs, please use the ctools Bug Tracker at
https://cta-redmine.irap.omp.eu/projects/ctools.  Before using the
tracker, please read
https://cta-redmine.irap.omp.eu/projects/ctools/wiki/Submission_guidelines


Known problems
==============

Solaris
-------

Although ctools compile on Solaris using the Sun compiler, there are
problems with global symbols in the shared GammaLib library that prevent
the model registry to work correctly.  Furthermore, GammaLib is not able
to catch its own exceptions, which prevents the FITS interface to work
correctly.  Possible it will work using gcc on Solaris, yet this has not
so far been tested.


Contact
=======
To get in touch with the ctools developers and to contribute to the
project please contact Juergen Knoedlseder <jurgen.knodlseder@irap.omp.eu>.
