.. _1dc_getting_data:

Getting the data
================

The data for the
:ref:`first CTA Data Challenge <glossary_1dc>`
are available on the
`CTA ownCloud server <https://owncloud.cta-observatory.org>`_.
You need your usual CTA credentials to download the data.

The data have been split in number of tarballs that are summarised in
the following table. You may click on the links in the second column
to download the data.

 +-------------------------------+-------------------------------------------------------------------------------+
 | Content                       | Download link                                                                 |
 +===============================+===============================================================================+
 | Instrument Response Functions | `<https://owncloud.cta-observatory.org/remote.php/webdav/1dc/caldb.tar.gz>`_  |
 +-------------------------------+-------------------------------------------------------------------------------+
 | Sky and background models     | `<https://owncloud.cta-observatory.org/remote.php/webdav/1dc/models.tar.gz>`_ |
 +-------------------------------+-------------------------------------------------------------------------------+

You may also download the files using the ``wget`` tool from the command
line by typing (``<user>`` needs to be replaced by your CTA user name):

.. code-block:: bash

   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/caldb.tar.gz --user=<user> --ask-password
   $ wget https://owncloud.cta-observatory.org/remote.php/webdav/1dc/models.tar.gz --user=<user> --ask-password

.. warning::
   Don't type the ``$`` symbol. It indicates in the following that a command
   should be typed on the command line of a console.

After downloading, uncompress the files at any place by typing

.. code-block:: bash

   $ tar xfvz caldb.tar.gz
   $ tar xfvz models.tar.gz

After that you should have the folder ``1dc`` in your current working
directory. Type there the following to set the ``CALDB`` and ``CTADATA``
environment variables:

.. code-block:: bash

   $ export CALDB=$PWD/1dc/caldb
   $ export CTADATA=$PWD/1dc/data

.. note::
   You may consider adding the ``CALDB`` and ``CTADATA`` environment variables
   to your ``.bashrc`` file (or equivalent) so that your analysis environment
   for the
   :ref:`first CTA Data Challenge <glossary_1dc>`
   is always setup correctly.

