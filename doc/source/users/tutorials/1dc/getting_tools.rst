.. _1dc_getting_tools:

Getting the analysis tools
==========================

We recommend that you use ctools to analyse the data from the
:ref:`first CTA Data Challenge <glossary_1dc>`.

The reference version for the pre-release of the data is ``ctools-1.2.0``.
``ctools-1.2.0`` depends on ``gammalib-1.2.0``. You may either download
the source code for both packages or use a pre-compiled binary package for
Mac OS X that includes all necessary code:

* `Mac OS X binary package <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.2.0-macosx10.7.dmg>`_ (Mac OS X 10.7+)
* `GammaLib <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.2.0.tar.gz>`_ source code tarball
* `ctools <http://cta.irap.omp.eu/ctools/releases/ctools/ctools-1.2.0.tar.gz>`_ source code tarball

After downloading, configure GammaLib and ctools as follows:

.. code-block:: bash

   $ export GAMMALIB=/usr/local/gamma
   $ source $GAMMALIB/bin/gammalib-init.sh
   $ export CTOOLS=/usr/local/gamma
   $ source $CTOOLS/bin/ctools-init.sh

Have a look at :ref:`getting the ctools <sec_getting>` if you need more
information about the installation and configuration process.

.. note::
   You may consider adding the ``GAMMALIB`` and ``CTOOLS`` environment variables
   as well as the sourcing of the configuration scripts to your ``.bashrc`` file
   (or equivalent) so that your analysis environment
   for the
   :ref:`first CTA Data Challenge <glossary_1dc>`
   is always setup correctly.



