.. _sec_install_devel:

Installing the development version
==================================

The current development version of the code can be downloaded as source
tarballs, Mac OS X binaries, or directly from the `Git <https://git-scm.com/>`_
repository.

The current ctools development release is ``ctools-1.5.0.dev1``.


Source tarballs
---------------

Download the source tarballs from the following links

- `ctools <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.5.0.dev1.tar.gz>`_
- `GammaLib <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.5.0.dev1.tar.gz>`_

and follow the instructions on :ref:`sec_install_source`.


Binary packages
---------------

Download the installer image from the following link

- `Mac OS X (10.7+) <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.5.0.dev1-macosx10.7.dmg>`_

and follow the instructions on :ref:`sec_install_binary`.


Git repository
--------------

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

before retrieving the code.
Before you will be able to compile the code you need to generate the
configuration file using the ``autogen.sh`` script.
Also make sure that you're actually on the devel branch of the git
repository. GammaLib and ctools can be compiled and configured using
the following command sequence (the code will be installed into the 
``/usr/local/gamma`` directory):

.. code-block:: bash

   $ cd gammalib
   $ git checkout devel
   $ ./autogen.sh
   $ ./configure
   $ make
   $ make check
   $ sudo make install
   $ export GAMMALIB=/usr/local/gamma
   $ source $GAMMALIB/bin/gammalib-init.sh
   $ cd ../ctools
   $ git checkout devel
   $ ./autogen.sh
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
   when fetching a release as the Python wrappers are bundled with the
   release tarballs.
