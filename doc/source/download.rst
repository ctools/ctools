.. _download:

Download
========

ctools can be obtained in form of releases or directly from the git 
development repository. Prefer a release if you intend using ctools
for production (and publications). Clone the code from git if you need
the most recent code that implements new features and corrects known
bugs.


Releases
--------

The latest ctools release is ``ctools-1.2.0`` (3 March 2017).

Below a table of ctools releases. Please note that at this stage of the
project there is a strict link between the ctools and gammalib versions.
Please make sure that you have the corresponding gammalib version installed
before installing ctools. The Mac OS X packages comprise both ctools and
gammalib.

.. list-table::
   :header-rows: 1
   :widths: 5 5 10

   * - `ctools <http://cta.irap.omp.eu/ctools/releases/ctools/ChangeLog>`_
     - `gammalib <http://cta.irap.omp.eu/ctools/releases/gammalib/ChangeLog>`_
     - Mac OS X package
   * - `1.2.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.2.0.tar.gz>`_
     - `1.2.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.2.0.tar.gz>`_
     - `ctools-1.2.0-macosx10.7.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.2.0-macosx10.7.dmg>`_
   * - `1.1.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.1.0.tar.gz>`_
     - `1.1.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.1.0.tar.gz>`_
     - `ctools-1.1.0-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.1.0-macosx10.3.dmg>`_
   * - `1.0.1 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.0.1.tar.gz>`_
     - `1.0.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.0.1.tar.gz>`_
     - `ctools-1.0.1-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.0.1-macosx10.3.dmg>`_
   * - `1.0.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.0.0.tar.gz>`_
     - `1.0.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.0.0.tar.gz>`_
     - `ctools-1.0.0-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.0.0-macosx10.3.dmg>`_
   * - `0.10.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-0.10.0.tar.gz>`_
     - `0.11.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-0.11.0.tar.gz>`_
     - `ctools-0.10.0-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-0.10.0-macosx10.3.dmg>`_
   * - `0.9.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-0.9.0.tar.gz>`_
     - `0.10.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-0.10.0.tar.gz>`_
     - `ctools-0.9.1-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-0.9.1-macosx10.3.dmg>`_
   * - `0.8.1 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-08-01.tar.gz>`_
     - `0.9.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-09-01.tar.gz>`_
     - `ctools-00-08-01-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-08-01-macosx10.3.dmg>`_
   * - `0.8.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-08-00.tar.gz>`_
     - `0.9.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-09-00.tar.gz>`_
     - `ctools-00-08-00-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-08-00-macosx10.3.dmg>`_
   * - `0.7.1 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-07-01.tar.gz>`_
     - `0.8.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-08-01.tar.gz>`_
     - `ctools-00-07-01-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-07-01-macosx10.3.dmg>`_
   * - `0.7.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-07-00.tar.gz>`_
     - `0.8.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-08-00.tar.gz>`_
     - `ctools-00-07-00-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-07-00-macosx10.3.dmg>`_
   * - `0.6.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-06-00.tar.gz>`_
     - `0.7.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-07-00.tar.gz>`_
     - `ctools-00-06-00-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-06-00-macosx10.3.dmg>`_
   * - `0.5.1 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-05-01.tar.gz>`_
     - `0.6.2 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-06-02.tar.gz>`_
     - `ctools-00-05-01-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-05-01-macosx10.3.dmg>`_
   * - `0.5.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-05-00.tar.gz>`_
     - `0.6.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-06-01.tar.gz>`_
     - `ctools-00-05-00-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-05-00-macosx10.3.dmg>`_
   * - `0.4.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-04-00.tar.gz>`_
     - `0.5.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-05-00.tar.gz>`_
     - `ctools-00-04-00-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-04-00-macosx10.3.dmg>`_
   * - `0.3.0 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-03-00.tar.gz>`_
     - `0.4.2 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-02.tar.gz>`_
     - `ctools-00-03-00-macosx10.3.dmg <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-03-00-macosx10.3.dmg>`_
   * - `0.2.5 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-02-05.tar.gz>`_
     - `0.4.11 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-11.tar.gz>`_
     -
   * - `0.2.4 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-02-04.tar.gz>`_
     - `0.4.10 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-10.tar.gz>`_
     -
   * - `0.2.3 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-02-03.tar.gz>`_
     - `0.4.9 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-09.tar.gz>`_
     -
   * - `0.2.1 <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-00-02-01.tar.gz>`_
     - `0.4.7 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-07.tar.gz>`_
     -


Development release
-------------------

The current ctools development release is ``ctools-1.3.0.dev1``.
This release reflects the status of the current ``devel`` branch of
the ctools git repository.

* `Mac OS X binary package <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.3.0.dev1-macosx10.7.dmg>`_ (Mac OS X 10.7+)
* `GammaLib <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.3.0.dev1.tar.gz>`_ source code tarball
* `ctools <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.3.0.dev1.tar.gz>`_ source code tarball


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

Please read the :ref:`sec_getting` section if you need more information on
how to install ctools.

.. note::

  You need `swig <http://www.swig.org/>`_ on your system to build the
  Python wrappers when you get the code from Git. Python wrappers are
  not stored in the Git repository but are built using
  `swig <http://www.swig.org/>`_ from interface files located in the
  pyext folder. However, you do not need `swig <http://www.swig.org/>`_
  when fetching a release as the Python wrappers are bundled with the
  release tarballs.
