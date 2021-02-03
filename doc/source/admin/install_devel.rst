.. _sec_install_devel:

Installing the development version
==================================

Alternatively to stable releases you can also install the development version
of the code which includes the latest features that will become available
with the next stable release. Installing the development version is recommended if
you want to use latest code features and if you cannot wait for the next
stable release. The current ctools development release is ``ctools-2.0.0.dev``.

Installing the development version is also necessary **if you want to contribute
to code development. In that case please follow the instructions on**
:ref:`using_git` **and ignore the rest of this page.**


Getting the development version
-------------------------------

The first step is to get the development version of the code. There are three
options for that first step: either get the source tarball or clone the
source code from the `Git <https://git-scm.com/>`_ repository or download
the binary installers for Mac OS.


Getting the source tarball
~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the source tarballs from the following links

- `ctools <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-2.0.0.dev.tar.gz>`_
- `GammaLib <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-2.0.0.dev.tar.gz>`_

Cloning the Git repository
~~~~~~~~~~~~~~~~~~~~~~~~~~

To clone the gammalib and ctools source codes, type

.. code-block:: bash

   $ git clone https://cta-gitlab.irap.omp.eu/gammalib/gammalib.git
   $ git clone https://cta-gitlab.irap.omp.eu/ctools/ctools.git
  
This will create directories named gammalib and ctools under the current
working directory that will contain the gammalib and ctools source code.
In case that the cloning does not work you may try adding

.. code-block:: bash

   $ export GIT_SSL_NO_VERIFY=true

or

.. code-block:: bash

   $ git config --global http.sslverify "false"

before retrieving the code. Before you will be able to compile the code directly
you need to generate the configuration file using the ``autogen.sh`` script.
Also make sure that you're actually on the devel branch of the git repository.
You do this by typing:

.. code-block:: bash

   $ cd gammalib
   $ git checkout devel
   $ ./autogen.sh
   $ cd ../ctools
   $ git checkout devel
   $ ./autogen.sh


Downloading Mac OS binary packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can download the Mac OS installer image from the following link

- `Mac OS (10.7+) <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-2.0.0.dev-macosx10.7.dmg>`_

and follow the instructions on :ref:`sec_install_binary`.


Compiling the development version
---------------------------------

The second step is to compile the development version. You do not need to
compile the code in case that you installed the Mac OS binary packages.

For code compilation, you can either compile the code directly using ``make``
or you can build a conda package. The latter is preferred if you want to use
ctools in a conda environment. We do not recommend to compile the development
version directly against a conda version of Python.


Compiling the development version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ctools and GammaLib can be compiled and configured using the following sequence
of commands (the code will be installed into the ``/usr/local/gamma`` directory):

.. code-block:: bash

   $ cd gammalib
   $ ./configure
   $ make
   $ make check
   $ sudo make install
   $ export GAMMALIB=/usr/local/gamma
   $ source $GAMMALIB/bin/gammalib-init.sh
   $ cd ../ctools
   $ ./configure
   $ make
   $ make check
   $ sudo make install
   $ export CTOOLS=/usr/local/gamma
   $ source $CTOOLS/bin/ctools-init.sh

Please read the :ref:`sec_install_source` section if you need more information
on how to install ctools.

.. note::
   You need `swig <http://www.swig.org/>`_ on your system to build the
   Python wrappers when you get the code from Git. Python wrappers are
   not stored in the Git repository but are built using
   `swig <http://www.swig.org/>`_ from interface files located in the
   pyext folder. However, you do not need `swig <http://www.swig.org/>`_
   when fetching a tarball as the Python wrappers are bundled with the
   tarballs.

Building a conda package
~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively you can create a conda package using the following sequence
of commands (make sure that anaconda is included in your ``$PATH`` environment):

.. code-block:: bash

   $ cd gammalib
   $ ./configure
   $ conda-build dev/conda.recipe
   $ cd ../ctools
   $ ./configure
   $ conda-build dev/conda.recipe

Once this is done, you can create a conda environment using the development
version as follows:

.. code-block:: bash

   $ conda create -n ctools-devel python=3.6
   $ source activate ctools-devel
   $ conda install --use-local ctools=2.0.0.dev

Note that you can choose between Python 2.7, 3.5, 3.6, 3.7, 3.8 and 3.9 for
your conda environment. If you only need the conda package for one specific
Python version you can build the conda package as follows:

.. code-block:: bash

   $ conda-build dev/conda.recipe --python 3.6

