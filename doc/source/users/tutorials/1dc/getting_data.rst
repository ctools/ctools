.. _1dc_getting_data:

Getting the data
================

The
:ref:`first CTA Data Challenge <glossary_1dc>`
data are available on the
`CTA ownCloud server <https://owncloud.cta-observatory.org>`_.
You need your usual CTA credentials to download the data.

The release contains simulated data for the Galactic Plane Survey, the
Galactic Centre Survey, the Extragalactic Survey and the AGN monitoring program,
:ref:`instrument response functions (IRFs) <glossary_irf>`, and the
:ref:`model definition files <glossary_moddef>`
that were used for the simulation.

The data were split into the separate files that can be downloaded
by clicking on the following links:

* `Galactic Plane Survey <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/gps.tar.gz>`_ (8.3 GB)
* `Galactic Centre Survey <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/gc.tar.gz>`_ (4.4 GB)
* `Extragalactic Survey <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/egal.tar.gz>`_ (2.5 GB)
* `AGN monitoring program <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/agn.tar.gz>`_ (??? GB)
* `Instrument Response Functions <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/caldb.tar.gz>`_ (1.2 MB)
* `Sky and background models <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/models.tar.gz>`_ (0.9 GB)

You can also download the files using the ``wget`` tool from the command
line by typing:

.. code-block:: bash

   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/gps.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/gc.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/egal.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/agn.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/caldb.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/models.tar.gz --user=<user> --ask-password

.. warning::
   ``<user>`` needs to be replaced by your CTA user name.
   Don't type the ``$`` symbol. It indicates that a command should be typed
   on the command line of your terminal.

After downloading, uncompress the files at any place by typing

.. code-block:: bash

   $ tar xfvz gps.tar.gz
   $ tar xfvz gc.tar.gz
   $ tar xfvz egal.tar.gz
   $ tar xfvz agn.tar.gz
   $ tar xfvz caldb.tar.gz
   $ tar xfvz models.tar.gz

You should now have a folder named ``1dc`` in your current working
directory with the following structure:

.. code-block:: bash

   1dc/
   1dc/caldb
   1dc/data
   1dc/models
   1dc/obs

Before continuing, please set the following environment variables:

.. code-block:: bash

   $ export CTADATA=$PWD/1dc
   $ export CALDB=$CTADATA/caldb

.. note::
   You may consider adding the ``CTADATA`` and ``CALDB`` environment variables
   to your ``.bashrc`` file (or equivalent) so that your analysis environment
   for the
   :ref:`first CTA Data Challenge <glossary_1dc>`
   is always setup correctly.

