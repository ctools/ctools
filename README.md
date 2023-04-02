ctools information
==================
* Version:             2.1.0.dev (2 April 2023)
* GammaLib dependency: 2.1.0.dev

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
ctools is a software package for the scientific analysis of astronomical
gamma-ray data. The software comprises an extensive set of tools for 
the analysis of data from existing and future Cherenkov telescopes, 
including H.E.S.S., VERITAS, MAGIC and CTA. ctools supports also the 
analysis of data from CGRO/COMPTEL, Fermi/LAT and INTEGRAL/SPI, 
enabling the exploration of the full gamma-ray energy band, spanning 
from hundreds of keV to hundreds of TeV.

The following tools and scripts are generic analysis utilities:

    ctbutterfly   - create a butterfly
    cterror       - likelihood profile error estimation
    ctmapcube     - generate sky map cube
    ctlike        - maximum likelihood model fitting
    cttsmap       - generate a TS map
    ctulimit      - compute upper limits
    cscaldb       - lists available instrument response functions
    csinfo        - checks ctools and GammaLib installations
    cslightcrv    - computes light curve
    csmodelinfo   - shows model container content
    csmodelmerge  - merges several model containers into one file
    csmodelselect - select models from model definition file
    csmodelsois   - generate map cube from subset of models
    csobsinfo     - shows observation container content
    csspec        - computes spectral points
    cstsdist      - generates Test Statistic distribution
    cstsmapmerge  - merges slices from Test Statistic map computations
    cstsmapsplit  - creates commands to split the Test Statistic map computations
    csworkflow    - run an analysis workflow

The following tools and script support the analysis of CTA and IACT data:

    ctbin         - event binning
    ctcubemask    - mask bins in binned analysis
    ctexpcube     - generate an exposure cube
    ctpsfcube     - generate a PSF cube
    ctedispcube   - generate energy dispersion cube
    ctbkgcube     - generate a background cube
    ctfindvar     - search for source variability
    ctmodel       - generation of model counts map
    ctobssim      - simulation of CTA observations
    ctphase       - computes the phase of each event
    ctprob        - computes event probability for a given model
    ctselect      - event selection
    ctskymap      - Sky mapping tool
    csbkgmodel    - generates background model for 3D analysis
    csebins       - generates energy boundaries for stacked analysis
    csobsdef      - generates observation definition file
    csobsselect   - select observations from observation definition file
    csphagen      - generates PHA, ARF, RMF files based on source/background regions
    csphasecrv    - computes phase curve
    cspull        - generates pull distribution
    csresmap      - generates residual map
    csresspec     - generates residual spectrum
    csroot2caldb  - creates a caldb entry from a ROOT file
    csscs         - Performs spectral component separation
    cssens        - computes CTA sensitivity
    cssrcdetect   - detects sources in sky map
    csviscube     - computes visibility cube

The following scripts support the management of an IACT database:

    csobs2caldb - Creates a caldb entry from an input observation
    csiactdata  - Shows information about IACT data available on the user machine
    csiactobs   - Generates observation definition file for IACT data from observation IDs
    csfindobs   - Generates a list of IACT observation IDs
    csiactcopy  - Copies IACT data from one location to another

The following scripts support COMPTEL science analysis:

    comlixfit    - Fit model to data using SRCLIX algorithm
    comlixmap    - Create TS map using SRCLIX algorithm
    comobsadd    - Combine observations
    comobsback   - Generate background model for COMPTEL observations
    comobsbin    - Bin COMPTEL observations
    comobsmodel  - Generate model for binned COMPTEL observations
    comobsres    - Generate residuals of COMPTEL observations
    comobsselect - Select observations from COMPTEL database
    comobssim    - Simulate COMPTEL observations
    comsrcdetect - Detect source in TS map


Web sites
=========
http://cta.irap.omp.eu/ctools                     - for ctools users
https://cta-redmine.irap.omp.eu/projects/ctools   - for ctools development


Prerequisites
=============
ctools require GammaLib.  Please refer to http://cta.irap.omp.eu/gammalib
for instructions about how to install GammaLib.

Once GammaLib is properly installed, make sure that you added the setup
script to your `.bashrc` or `$HOME/.profile` script:

    export GAMMALIB=/usr/local/gamma
    source $GAMMALIB/bin/gammalib-init.sh

If you use C shell or a variant then add the following to your
.cshrc or .tcshrc script:

    setenv GAMMALIB /usr/local/gamma
    source $GAMMALIB/bin/gammalib-init.csh

If you have installed GammaLib in another directory than `/usr/local/gamma`,
please adapt the path correspondingly.

If you really insist, you may install ctools in a directory different to
that hosting GammaLib, but we highly recommend to install both packages
together.


Conda Installation
==================
The easiest is to install ctools via conda.  This also takes care of the
installation of GammaLib.  Assuming that you have installed anaconda, type
the following:

    $ conda config --append channels conda-forge
    $ conda config --append channels cta-observatory
    $ conda install ctools


Linux Installation
==================
To build and install ctools, simply type the following:

    $ ./configure
    $ make
    $ make check
    $ make install

If the folder does not contain any `configure` file, please run

    $ ./autogen.sh 

before invoking `configure`.

By default ctools installs itself in `/usr/local/gamma`.  If you need to
install ctools in a different location or in your home directory, use
the `--prefix` option to `./configure`.  For example:

    $ ./configure --prefix=/home/yourname/projects
    $ make
    $ make check
    $ make install

The file INSTALL details more about using configure. Also try

    $ ./configure --help.

The `make check` command will run an extensive unit test to verify that
ctools was correctly built.  Make sure that all tests were successful. 


Macintosh OS Installation
==========================
ctools builds and installs seamlessly on all Mac OS starting from at least Mac
OS 10.6.

On older systems you may use Mac specific configure options:

    $ ./configure --enable-universalsdk[=PATH]

creates a universal build of ctools.  The optional argument specifies
which MacOS SDK should be used to perform the build.  This defaults to
`/Developer/SDKs/MacOSX.10.4u.sdk`.  Specify `/` when building on a 10.5
system or higher, especially when building 64-bit code.

    $ ./configure --with-univeral-archs=VALUE

specifies the kind of universal build that should be created.  Possible
values are: `32-bit`, `3-way`, `intel` or `all`.  By default, a `32-bit` 
build will be made.  This option is only valid when
`--enable-universalsdk` is specified.

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
ctools has been successfully installed on FreeBSD.  Make sure that you
use `gmake` for building of the software.  Otherwise follow the Linux
installation instructions above.


Solaris Installation
====================
ctools compile on Solairs, but there is an issue with using it as a
shared library (see "Known problems" below).


Windows Installation
====================
On Windows ctools needs to be installed into a virtual machine running
a Linux distribution.


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
been tested so far.


Contact
=======
To get in touch with the ctools developers and to contribute to the
project please contact Juergen Knoedlseder <jurgen.knodlseder@irap.omp.eu>.
