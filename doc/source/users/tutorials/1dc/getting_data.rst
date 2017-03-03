.. _1dc_getting_data:

Getting the data
================

A pre-release of the data for the
:ref:`first CTA Data Challenge <glossary_1dc>`
is available on the
`CTA ownCloud server <https://owncloud.cta-observatory.org>`_.
You need your usual CTA credentials to download the data.
The full
:ref:`first CTA Data Challenge <glossary_1dc>`
data release is planned for the second half of March 2017.

The pre-release contains simulated data for the Galactic Centre Survey,
:ref:`Instrument Response Functions (IRFs) <glossary_irf>`, and the
:ref:`model definition files <glossary_moddef>`
that were used for the simulation.
The data were split into the separate files that can be downloaded
by clicking on the following links:

* `Galactic Centre Survey <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.pre/gc-pre.tar.gz>`_ (5.9 GB)
* `Instrument Response Functions <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.pre/caldb-pre.tar.gz>`_ (0.2 MB)
* `Sky and background models <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.pre/models-pre.tar.gz>`_ (0.9 GB)

You can also download the files using the ``wget`` tool from the command
line by typing:

.. code-block:: bash

   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.pre/gc-pre.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.pre/caldb-pre.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.pre/models-pre.tar.gz --user=<user> --ask-password

.. warning::
   ``<user>`` needs to be replaced by your CTA user name.
   Don't type the ``$`` symbol. It indicates that a command should be typed
   on the command line of a console.

After downloading, uncompress the files at any place by typing

.. code-block:: bash

   $ tar xfvz gc-pre.tar.gz
   $ tar xfvz caldb-pre.tar.gz
   $ tar xfvz models-pre.tar.gz

You should now have a folder named ``1dc.pre`` in your current working
directory with the following structure:

.. code-block:: bash

   1dc.pre/
   1dc.pre/caldb
   1dc.pre/data
   1dc.pre/models
   1dc.pre/obs

Before continuing, please set the following environment variables:

.. code-block:: bash

   $ export CTADATA=$PWD/1dc.pre
   $ export CALDB=$CTADATA/caldb

.. note::
   You may consider adding the ``CTADATA`` and ``CALDB`` environment variables
   to your ``.bashrc`` file (or equivalent) so that your analysis environment
   for the
   :ref:`first CTA Data Challenge <glossary_1dc>`
   is always setup correctly.

