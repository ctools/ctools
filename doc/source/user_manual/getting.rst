Getting the ctools
==================

Before you start
----------------

The procedure for building and installing ctools is modeled on GNU
software distributions. You will need the following to build the
software:

-  About 25 MB of free disk space.

-  An ANSI C++ compiler. ctools builds well using GNU g++ or clang.

-  make, automake, autoconf, libtools

-  gammalib (see next section).

-  Python, including the Python developer package that provides the 
   ``Python.h`` header files. Although ctools compiles without Python and/or
   the Python developer package installed, compiling with Python support
   is highly recommended to enable dynamic scripting of ctools and
   developing and using of cscripts.

No system administrator privileges are needed to install ctools.


Installing gammalib
-------------------

ctools is built on top of gammalib, hence gammalib needs to be built,
installed and configured before ctools can be installed. In case that
ctools is already installed on your system you may skip reading this
section and continue with Section TBD that explains how to install ctools.

You will need the following to build the gammalib:

-  About 200 MB of free disk space.

-  An ANSI C++ compiler. We recommend building GammaLib with the GNU g++
   or the clang compilers.

-  make, automake, autoconf, libtools

-  The cfitsio library for FITS file support together with the developer
   package that provides the ``cfitsio.h`` header.

-  Python, including the Python developer package that provides the 
   ``Python.h`` header files. Although ctools compiles without Python and/or
   the Python developer package installed, compiling with Python support
   is highly recommended to enable dynamic scripting of ctools and
   developing and using of cscripts.

Furthermore it is recommended to have readline and ncurses installed
(including the developer packages).

After downloading the gammalib tarball (see :ref:`download`), save it
in an appropriate location (for example ``$HOME/builds``), and type

.. code-block:: bash

  $ tar xvfz gammalib-1.0.0.tar.gz

(the ``$`` symbol indicates the console prompt and is not part of the
command that you should type in).

Step in the created directory and build, check and install gammalib by
typing

.. code-block:: bash

  $ cd gammalib-1.0.0
  $ ./configure
  $ make
  $ make check
  $ make install

The last step may require system administator privileges. In this case,
use

  $ sudo make install

If you do not have system administrator privileges, change the install
directory by typing for example

.. code-block:: bash

  $ ./configure --prefix=$HOME/gamma

which then installs GammaLib in the ``gamma`` folder of your home
directory.


Setting up the GammaLib environment
-----------------------------------

Before building and installing ctools, you have to configure GammaLib by
setting up some environment variables. This will be done automatically by
an initialisation script that is found in the ``bin`` directory of the 
GammaLib directory (we call here "GammaLib directory" the directory into
which GammaLib has been installed). Assuming that you have installed 
GammaLib in the default directory ``/usr/local/gamma`` you need to add the
following to your ``$HOME/.bashrc`` or ``$HOME/.profile`` script on a Linux 
machine:

.. code-block:: bash

  export GAMMALIB=/usr/local/gamma
  source $GAMMALIB/bin/gammalib-init.sh

If you use C shell or a variant then add the following to your 
``$HOME/.cshrc`` or ``$HOME/.tcshrc`` script:

.. code-block:: csh

  setenv GAMMALIB /usr/local/gamma
  source $GAMMALIB/bin/gammalib-init.csh


Building ctools
---------------

After downloading the ctools tarball (see :ref:`download`), save it in 
an appropriate location (for example ``$HOME/builds``), and type

.. code-block:: bash

  $ tar xvfz ctools-1.0.0.tar.gz

(the ``$`` symbol indicates the console prompt and is not part of the
command that you should type in).

Step in the directory and build the ctools by typing

.. code-block:: bash

  $ cd ctools-1.0.0
  $ ./configure
  $ make

at the operating system prompt. The ``./configure`` command customizes
the Makefiles for the particular system, the ``make`` command compiles
the source files and builds the executables. Type ``./configure`` and
not simply ``configure`` to ensure that the configuration script in the
current directory is run and not some other system-wide configuration script. 

You can get the full list of configuration options by typing

.. code-block:: bash

  $ ./configure --help


Testing ctools
--------------

Before installing the ctools you should execute the unit test suite to 
make sure that ctools have been built correctly. For this, type

.. code-block:: bash

  $ make check

If you have automake version 1.13 or newer installed, you should see the
following output at the end of the unit testing:

.. code-block:: bash

  PASS: test_ctools.sh
  PASS: test_cscripts.sh
  PASS: test_python.py
  make[4]: Nothing to be done for `all'.
  ============================================================================
  Testsuite summary for ctools 1.0.0
  ============================================================================
  # TOTAL: 3
  # PASS:  3
  # SKIP:  0
  # XFAIL: 0
  # FAIL:  0
  # XPASS: 0
  # ERROR: 0
  ============================================================================

For older automake version, you should see

.. code-block:: bash

  ***************
  * Test ctools *
  ***************
  Test ctobssim: .... ok
  Test ctskymap: .. ok
  Test ctbin: .... ok
  Test ctmodel: ...... ok
  Test ctselect: .. ok
  Test ctlike: ........ ok
  Test cttsmap: .. ok
  Test ctexpcube: .... ok
  Test ctpsfcube: .... ok
  Test ctbkgcube: ...... ok
  Test ctcubemask: .. ok
  Test ctbutterfly: .. ok
  PASS: test_ctools.sh

  *****************
  * Test cscripts *
  *****************
  Test cspull: .. ok
  Test cstsdist: .. ok
  Test csresmap: .. ok
  PASS: test_cscripts.sh

  ***********************
  * ctools unit testing *
  ***********************
  Test ctobssim functionality: .......... ok
  Test ctobssim on observation container: .... ok
  Test ctselect functionality: ... ok
  Test ctbin functionality: ... ok
  Test ctlike functionality: ... ok
  Test cttsmap functionality: ... ok
  Test ctmodel functionality: ... ok
  Test ctskymap functionality: ... ok
  Test ctexpcube functionality: ... ok
  Test ctpsfcube functionality: ... ok
  Test ctbkgcube functionality: ... ok
  Test ctcubemask functionality: ... ok
  Test ctbutterfly functionality: ... ok
  Test unbinned pipeline with FITS file saving: .... ok
  Test unbinned in-memory pipeline: .... ok
  PASS: test_python.py
  ==================
  All 3 tests passed
  ==================

The same detailed information is also available for the newer automake 
versions, but there it is written in log files that you can find in the 
``test`` directory of the ctools:

.. code-block:: bash

  test_ctools.sh.log
  test_cscripts.sh.log
  test_python.py.log

If you do not see the same output, but a failure message, please report 
this to the ctools developer team.


Installing ctools
-----------------

The ctools are installed by typing

.. code-block:: bash

  $ make install

If the destination directory is owned by ``root`` (which is normally the 
case when using the default), administrator privileges are needed for
installation. In this case, type

.. code-block:: bash

  $ sudo make install

By default, the install directory is set to ``/usr/local/gamma``. To 
change the install directory (for example in case that you do not
have system administrator privileges), an optional ``--prefix`` argument
can be given, for example:

.. code-block:: bash

  $ ./configure --prefix=$HOME/gamma


Setting up the ctools environment
---------------------------------

You have to configure ctools by setting up some environment variables. This
will be done automatically by an initialisation script that is found in the
``bin`` directory of the directory into which ctools has been installed. 
Assuming that you have installed ctools into ``/usr/local/gamma`` you need
to add the following to your ``$HOME/.bashrc`` or ``$HOME/.profile`` script
on a Linux machine:

.. code-block:: bash

  export CTOOLS=/usr/local/gamma
  source $CTOOLS/bin/ctools-init.sh

If you use C shell or a variant then add the following to your 
``$HOME/.cshrc`` or ``$HOME/.tcshrc`` script:

.. code-block:: csh

  setenv CTOOLS /usr/local/gamma
  source $CTOOLS/bin/ctools-init.csh


Known problems
--------------

GammaLib unit tests fail
~~~~~~~~~~~~~~~~~~~~~~~~

Some users have reported failure of a large fraction of the GammaLib unit
tests after after typing ``make check``. In all cases, this was related to
the absence of the directory where the shared ``libcfitsio`` library 
resides in the library load path. To solve the issue, locate the directory
where the shared ``libcfitsio`` library resides and then type

.. code-block:: bash

  export LD_LIBRARY_PATH=/directory/to/lib:$LD_LIBRARY_PATH

on Unix based systems or

.. code-block:: bash

  export DYLD_LIBRARY_PATH=/directory/to/lib:$DYLD_LIBRARY_PATH

on Mac OS X (``/directory/to/lib`` should be replaced by the correct
library path on your system).


Python support
~~~~~~~~~~~~~~

ctools comes with Python wrappers so that all classes can be directly 
used from Python. To compile-in Python support, ctools needs the 
``Python.h`` header file, which on many distributions is not installed
by default. To make ``Python.h`` available, install the Python developer
package in your distribution. Otherwise you will not be able to use ctools
from Python.

GammaLib and ctools release bundles come with Python wrapper files, but 
these wrapper files will be deleted once you type ``make clean``. Wrapper
files will be built automatically when needed, but this required
`swig <http://www.swig.org/>`_ installed on your system. In that case,
please make sure that you swig version 2.0.4 or newer.


Mac OS X
~~~~~~~~

The Python development package is not installed by default on Mac OS X,
and consequently, the ``Python.h`` header file is missing that is needed
to compile in the Python bindings. The configure script recognises this
fact and adjust the build procedure accordingly, but you will not be able
to use ctools as a Python module. So better install the Python development
package before installing ctools.

It was also reported that adding of 
``export MACOSX_DEPLOYMENT_TARGET=10.6`` to the ``.bashrc`` file was
necessary one some Mac OS X 10.6 installations to make ctools compile.


Solaris
~~~~~~~

Although ctools builds on Solaris using the Sun compiler, there are problems
with global symbols in shared libraries and exception catching, which prevents
the FITS interface to work correctly. ctools has however been built and tested
successfully using the GNU compiler, and this is the only build method that
is currently supported. Problems have also been encountered when compiling
cfitsio versions more recent than 3.250. The problems have been reported to
the cfitsio developer team, and are likely to be solved in the future. For 
the time being, it is recommended to use cfitsio version 3.250 on Solaris.


OpenSolaris
~~~~~~~~~~~

On OpenSolaris, the same problems concerning the SunStudio compiler occur
as for Solaris, and also here, the GNU compiler is the recommended tool to
build ctools. Also here, cfitsio version 3.250 is the recommended library as
more recent version feature relocation problems. ctools has been tested 
using gcc 4.3.2 on OpenSolaris 2009.06. Make sure to create the symbolic 
links

.. code-block:: csh

  $ ln -s /usr/bin/gcc4.3.2 /usr/bin/gcc
  $ ln -s /usr/bin/g++4.3.2 /usr/bin/g++

which are not there by default to avoid excess warnings during compilation.