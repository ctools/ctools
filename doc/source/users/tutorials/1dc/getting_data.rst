.. _1dc_getting_data:

Getting the data
================

The
:ref:`first CTA Data Challenge <glossary_1dc>`
data are available on the
`CTA ownCloud server <https://owncloud.cta-observatory.org>`_.
Only CTA Consortium members are entitled to access the data.
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
* `AGN monitoring program <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/agn.wobble.tar.gz>`_ (4.7 GB)
* `Instrument Response Functions <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/caldb.tar.gz>`_ (1.2 MB)
* `Sky and background models <https://owncloud.cta-observatory.org/remote.php/webdav/1dc/models.tar.gz>`_ (0.9 GB)

You can also download the files using the ``wget`` tool from the command
line by typing:

.. code-block:: bash

   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/gps.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/gc.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/egal.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/agn.wobble.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/caldb.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/models.tar.gz --user=<user> --ask-password

.. warning::
   ``<user>`` needs to be replaced by your CTA user name.
   Don't type the ``$`` symbol. It indicates that a command should be typed
   on the command line of your terminal.

.. note::
   You may need to add the ``--no-check-certificate`` option to download the
   files using ``wget``.

.. note::
   In case you encounter the following error message

   .. code-block:: bash

      HTTP request sent, awaiting response... 401 Unauthorized
      Authentication selected: Basic realm="ownCloud"
      Reusing existing connection to owncloud.cta-observatory.org:443.
      HTTP request sent, awaiting response... 404 Not Found
      2017-06-06 14:43:33 ERROR 404: Not Found.

   or any other problems with the download, please contact
   sstober@cta-observatory.org.

To make sure that you have the relevant files you may check the MD5 checksums
after downloading by typing for example

.. code-block:: bash

   $ md5 gps.tar.gz

on Mac OS X or

.. code-block:: bash

   $ md5sum gps.tar.gz

on Linux. Below the expected results:

.. code-block:: bash

   gps.tar.gz         c7129c71e911920f82273270c9d7c818
   gc.tar.gz          07aa46bb2d67cd64fa4685078643b043
   egal.tar.gz        4ea5f2ff8312b700dcc133a89d9fd096
   agn.wobble.tar.gz  df8fc4f9de22ac36a86023ccbfd0cf89
   caldb.tar.gz       ce718128e1b2fe19c4b8dc79b688d113
   models.tar.gz      87662cf53a2988522f485cad984d8810

Now uncompress the files at any place by typing

.. code-block:: bash

   $ tar xfvz gps.tar.gz
   $ tar xfvz gc.tar.gz
   $ tar xfvz egal.tar.gz
   $ tar xfvz agn.wobble.tar.gz
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

