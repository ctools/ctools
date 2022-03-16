.. _1dc_getting_data:

Getting the data
================

The data for the
:ref:`first CTA Data Challenge <glossary_1dc>`
are only available to members of the CTA Consortium.
Instructions for CTA Consortium of how to download the data can be found
`here <https://forge.in2p3.fr/projects/data-challenge-1-dc-1/wiki/Getting_data>`_
(access to the instructions are restricted).

Once downloaded, uncompress the files at any place by typing

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

