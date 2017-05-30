.. _1dc_getting_data:

Getting the data
================

The
:ref:`first CTA Data Challenge <glossary_1dc>`
data are available on the
`CTA ownCloud server <https://owncloud.cta-observatory.org>`_.
You need your usual CTA credentials to download the data.

The release contains simulated data for the Galactic Plane Survey, the
Galactic Centre Survey and the Extragalactic Survey,
:ref:`Instrument Response Functions (IRFs) <glossary_irf>`, and the
:ref:`model definition files <glossary_moddef>`
that were used for the simulation.
In addition, ctools metadata, including
:ref:`Observation Definition Files <glossary_obsdef>`
are available.

The data were split into the separate files that can be downloaded
by clicking on the following links:

* `Galactic Plane Survey <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/gps-south.tar.gz>`_ (6.1 GB)
* `Galactic Centre Survey <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/gc-south.tar.gz>`_ (4.9 GB)
* `Extragalactic Survey <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/egal-south.tar.gz>`_ (1.1 GB)
* `Instrument Response Functions <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/caldb-south.tar.gz>`_ (0.1 MB)
* `Sky and background models <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/models-south.tar.gz>`_ (0.9 GB)
* `ctools metadata <https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/obs-south.tar.gz>`_ (0.2 MB)

You can also download the files using the ``wget`` tool from the command
line by typing:

.. code-block:: bash

   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/gps-south.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/gc-south.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/egal-south.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/caldb-south.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/models-south.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc.south/obs-south.tar.gz --user=<user> --ask-password

.. warning::
   ``<user>`` needs to be replaced by your CTA user name.
   Don't type the ``$`` symbol. It indicates that a command should be typed
   on the command line of your terminal.

After downloading, uncompress the files at any place by typing

.. code-block:: bash

   $ tar xfvz gps-south.tar.gz
   $ tar xfvz gc-south.tar.gz
   $ tar xfvz egal-south.tar.gz
   $ tar xfvz caldb-south.tar.gz
   $ tar xfvz models-south.tar.gz
   $ tar xfvz obs-south.tar.gz

You should now have a folder named ``1dc.south`` in your current working
directory with the following structure:

.. code-block:: bash

   1dc.south/
   1dc.south/caldb
   1dc.south/data
   1dc.south/models
   1dc.south/obs

Before continuing, please set the following environment variables:

.. code-block:: bash

   $ export CTADATA=$PWD/1dc.south
   $ export CALDB=$CTADATA/caldb

.. note::
   You may consider adding the ``CTADATA`` and ``CALDB`` environment variables
   to your ``.bashrc`` file (or equivalent) so that your analysis environment
   for the
   :ref:`first CTA Data Challenge <glossary_1dc>`
   is always setup correctly.

