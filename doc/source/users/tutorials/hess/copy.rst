.. _iact_copy:

Copying IACT data
==================

  .. admonition:: What you will learn

     You will learn how copy data from an IACT data store to your disk.

Before you start
----------------

Currently the way how IACT data is stored is quite prototypical.
The ctools can work with IACT data that is stored following the specifications
provided in the
`Data formats for gamma-ray astronomy <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/index.html>`_.
The script :ref:`csiactcopy` can be used to transfer IACT data from one system
to another.
IACT data is the property of collaborations (such as HESS, MAGIC or VERITAS)
and not publicly available, except for the recent
`public release of some H.E.S.S. data <https://www.mpi-hd.mpg.de/hfm/HESS/pages/dl3-dr1>`_.
In this tutorial we use these public data.
After downloading the data and uncompression the tarball you should have a
folder named ``hess_dl3_dr1`` in your current workingcdirectory with the
following structure:

.. code-block:: bash

   hess_dl3_dr1/
   hess_dl3_dr1/data
   hess_dl3_dr1/hdu-index.fits.gz
   hess_dl3_dr1/obs-index.fits.gz
   hess_dl3_dr1/README.txt

You have to add the file ``master.json`` to this folder with the following
content:

.. code-block:: json

   {
     "datasets": [
       {
         "name": "hess_dl3_dr1",
         "hduindx": "hdu-index.fits.gz",
         "obsindx": "obs-index.fits.gz",
         "Total observations": 106,
         "Total time": "48.6 hours",
         "Total events": 662826,
         "Disk usage": "43 MB",
         "Date created": "2018-09-10"
       }
     ]
   }

.. note::
   If you are member of an IACT collaboration and if the collaboration has
   setup a datastore in the
   `compliant format <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/index.html>`_
   you can also establish a remote connection to this datastore using for
   example

   .. code-block:: bash

      $ sshfs user@remote.server.com:/path/to/remote/fits/data/ /project-data/hess

   The success of this operation can be verified by checking if there is a file
   called ``master.json`` in the mounted directory.


Check what data are available
-----------------------------

Before you start an analysis, it is important that you know what kind of data
are available for analysis.
In the `current scheme <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/index.html>`_
every dataset is represented by a unique name, in our case ``hess_dl3_dr1``.
You have to pass this string to some set of scripts to find and use the data
for analysis.

To check what datasets are available you should run the :ref:`csiactdata`
script:

.. code-block:: bash

   $ csiactdata debug=yes
   Path were data is located [] /project-data/hess/hess_dl3_dr1

This script will list the names of available datasets. Here at excerpt of the
log file ``csiactdata.log``:

.. code-block:: none

   2019-04-05T13:10:35: +======================+
   2019-04-05T13:10:35: | Data storage entries |
   2019-04-05T13:10:35: +======================+
   2019-04-05T13:10:35:  Master index file .........: /project-data/hess/hess_dl3_dr1/master.json
   2019-04-05T13:10:35:
   2019-04-05T13:10:35: +------------------------+
   2019-04-05T13:10:35: | Available data configs |
   2019-04-05T13:10:35: +------------------------+
   2019-04-05T13:10:35: === hess_dl3_dr1 ===
   2019-04-05T13:10:35:  Name ......................: hess_dl3_dr1
   2019-04-05T13:10:35:  Observation index .........: obs-index.fits.gz
   2019-04-05T13:10:35:  HDU index .................: hdu-index.fits.gz
   2019-04-05T13:10:35:  Date created ..............: 2018-09-10
   2019-04-05T13:10:35:  Total events ..............: 662826
   2019-04-05T13:10:35:  Total time ................: 48.6 hours
   2019-04-05T13:10:35:  Disk usage ................: 43 MB
   2019-04-05T13:10:35:  Total observations ........: 105


Copy data 
---------

You can now copy the data from the data store to your disk. In this example we
assume that the H.E.S.S. datastore is at ``/project-data/hess/hess_dl3_dr1``.
To do the copy you should type

.. code-block:: bash

   $ csiactcopy
   Location of remote master file [/path/to/mountpoint/master.json] /project-data/hess/hess_dl3_dr1/master.json
   Name of FITS production to download [iact-fits] hess_dl3_dr1
   Destination path of FITS data [/path/to/local/fits/data/] data

.. note::
   To monitor the progress on the screen you should set the hidden parameter
   ``debug=yes``. You may also increase the chattiness of the output by
   increasing the hidden parameter ``chatter`` from its default value of 2
   to for example 3.

.. note::
   If you are only interested in a subset of the dataset you can also download
   only that subset. For that purpose you have to specify an ASCII file
   containing a list of observation identifiers to be downloaded, e.g.

   .. code-block:: none

      20136
      20137
      20151

   You can then run :ref:`csiactcopy` by specifying the ASCII filename
   for the hidden ``runlist`` parameter:

   .. code-block:: bash

      $ csiactcopy runlist=runlist.txt
      Location of remote master file [/project-data/hess/hess_dl3_dr1/master.json]
      Name of FITS production to download [hess_dl3_dr1]
      Destination path of FITS data [data]


Troubleshooting
---------------
In case a download was interupted (e.g. the remote file system was
disconnected), the script can simply be rerun.
The index files, which list the files available on the system, are refreshed
after the copying process.
The hidden boolean parameter ``clobber`` (by default ``yes``) specifies if the
user wishes that local content should be overwritten with remote content.
In turn, specifying ``clobber=no`` can speed up the process since remote
content does not need to be copied if it is already available on the user
machine.
