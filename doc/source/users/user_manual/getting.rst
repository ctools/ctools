.. _sec_getting:

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


.. _sec_installing_gammalib:

Installing GammaLib
-------------------

ctools is built on top of GammaLib, hence GammaLib needs to be built,
installed and configured before ctools can be installed. In case that
ctools is already installed on your system you may skip reading this
section and continue with the :ref:`sec_building_ctools` section.

You will need the following to build the GammaLib:

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

After downloading the GammaLib tarball (see :ref:`download`), save it
in an appropriate location (for example ``$HOME/builds``), and type

.. code-block:: bash

  $ tar xvfz gammalib-1.2.0.tar.gz

(the ``$`` symbol indicates the console prompt and is not part of the
command that you should type in).

Step in the created directory and build, check and install gammalib by
typing

.. code-block:: bash

  $ cd gammalib-1.2.0
  $ ./configure
  $ make
  $ make check
  $ make install

The last step may require system administator privileges. In this case,
use

.. code-block:: bash

  $ sudo make install

If you do not have system administrator privileges, change the install
directory by typing for example

.. code-block:: bash

  $ ./configure --prefix=$HOME/gamma

which then installs GammaLib in the ``gamma`` folder of your home
directory.


.. _sec_setup_gammalib:

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


.. _sec_building_ctools:

Building ctools
---------------

After downloading the ctools tarball (see :ref:`download`), save it in 
an appropriate location (for example ``$HOME/builds``), and type

.. code-block:: bash

  $ tar xvfz ctools-1.2.0.tar.gz

(the ``$`` symbol indicates the console prompt and is not part of the
command that you should type in).

Step in the directory and build the ctools by typing

.. code-block:: bash

  $ cd ctools-1.2.0
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


.. _sec_testing_ctools:

Testing ctools
--------------

Before installing the ctools you should execute the unit test suite to 
make sure that ctools have been built correctly. For this, type

.. code-block:: bash

  $ make check

If you have automake version 1.13 or newer installed, you should see the
following output at the end of the unit testing:

.. code-block:: bash

   PASS: test_python_ctools.sh
   PASS: test_python_cscripts.sh
   PASS: test_examples.py
   ============================================================================
   Testsuite summary for ctools 1.2.0
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

   ***********************
   * ctools unit testing *
   ***********************
   Test ctobssim on command line: ..... ok
   Test ctobssim from Python: ................................................. ok
   Test ctselect on command line: ....... ok
   Test ctselect from Python: ..... ok
   Test ctbin on command line: ..... ok
   Test ctbin from Python: ............................................. ok
   ...
   PASS test_python_ctools.sh

   *************************
   * cscripts unit testing *
   *************************
   Test cscaldb on command line: .. ok
   Test cscaldb from Python: . ok
   Test cslightcrv on command line: ............... ok
   Test cslightcrv from Python: ................................................. ok
   Test csmodelinfo on command line: ..... ok
   Test csmodelinfo from Python: ..... ok
   ...
   PASS test_python_cscripts.sh

   ********************
   * Examples testing *
   ********************
   Test generate_prod3_irfs.py: .. ok
   Test make_pointings.py: ....... ok
   Test make_spectrum.py: .. ok
   Test make_ts_distributions.py: .. ok
   ...
   PASS test_examples.py
   ==================
   All 3 tests passed
   ==================

The same detailed information is also available for the newer automake 
versions, but there it is written in log files that you can find in the 
``test`` directory of the ctools:

.. code-block:: bash

  test_python_ctools.sh.log
  test_python_cscripts.sh.log
  test_examples.py.log

If you do not see the same output, but a failure message, please check
first the :ref:`issues` section. If you cannot fix the problem, please
create an issue on the ctools tracker
`here <https://cta-redmine.irap.omp.eu/projects/ctools>`_.

.. _sec_installing_ctools:

Installing ctools
-----------------

Now you are ready to install the ctools by typing

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

.. _sec_setup_ctools:

Setting up the ctools environment
---------------------------------

You have to configure ctools by setting up some environment variables. This
will be done automatically by an initialisation script that is found in the
``bin`` directory of the ctools installation. 
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

.. _sec_check_ctools_setup:

Checking your setup
-------------------

Now you should be ready to get started using Gammalib and ctools.

As a quick check that your setup is okay you can run ``csinfo check``:

.. code-block:: bash

  $ csinfo check

  Gammalib / ctools setup check:

     GAMMALIB environment variable ... ok
     CTOOLS   environment variable ... ok
     gammalib Python import .......... ok
     ctools   Python import .......... ok
     cscripts Python import .......... ok

     ===> Your Gammalib / ctools setup is OK.

If the setup is not okay, run the ``csinfo info`` command to print
detailed information about your setup. There's also a ``csinfo list``
command to quickly list the available tools.


.. _sec_known_problems:

Known problems
--------------

In case you encounter problem, please check the list of known
:ref:`installation_issues`.
If you encounter problems during the GammaLib installation, please
check the list of
`known GammaLib issues <http://cta.irap.omp.eu/gammalib/doc/html/issues.html>`_.
If you cannot solve your problems, please create an issue on the
ctools tracker
`here <https://cta-redmine.irap.omp.eu/projects/ctools>`_.
