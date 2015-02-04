Download
========

ctools can be obtained in form of releases or directly from the git 
development repository. Prefer a release if you intend using ctools
for production (and publications). Clone the code from git if you need
the most recent code that implements new features and corrects known
bugs.


Releases
--------

The latest ctools release is ``ctools-00-08-01`` (8 January 2015).

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

Git repository
--------------

To clone the ctools source code, type

.. code-block:: bash

  $ git clone https://cta-git.irap.omp.eu/ctools
  
This will create a directory named ctools under the current working directory
that will contain the ctools source code. Before you are able to compile the
code you need to generate the configuration file (make sure that you're
actually on the devel branch of the git repository):

.. code-block:: bash

  $ cd ctools
  $ git checkout devel
  $ ./autogen.sh

.. note::

  You need `swig <http://www.swig.org/>`_ on your system to build the
  Python wrappers when you get the code from Git. Python wrappers are
  not stored in the Git repository but are built using
  `swig <http://www.swig.org/>`_ from interface files located in the
  pyext folder. However, you do not need `swig <http://www.swig.org/>`_
  when fetching a release as the Python wrappers are bundled with the
  release tarballs.
