.. _1dc_getting_data:

Getting the data
================

The data for the
:ref:`First CTA Data Challenge <glossary_1dc>`
are bundled with the
:ref:`Instrument Response Functions <glossary_irf>`
into a single tarball.
Click on
`1dc.tar.gz <https://www.dropbox.com/s/2aru9gak10h9su6/1dc.tar.gz?dl=0>`_
to download the tarball (file size: 485 MB).
You can then uncompress and deploy the data at any place by typing

.. code-block:: bash

   $ tar xfvz 1dc.tar.gz

.. warning::
   Don't type the ``$`` symbol. It indicates in the following that a command
   should be typed on the command line of a console.

After deploying you should have the folder ``1dc`` in your current working
directory. Within that directory, please set the following environment variables:

.. code-block:: bash

   $ export CALDB=$PWD/1dc/caldb
   $ export CTADATA=$PWD/1dc/data
