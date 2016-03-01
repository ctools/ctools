.. _sec_copy:

Copying IACT data
==================

Before you start
----------------
Currently the way how IACT data is stored is quite prototypical. The ctools can work
with IACT data that is stored following the `specifications here <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/index.html>`_.
The script :ref:`csiactcopy` can be used to transfer these data from one system to another.
IACT data is the property of collaborations (such as HESS, MAGIC or VERITAS) and not publicly available. If you are member of one of these
collaborations, you should have access to a central place where the FITS data is stored and maintained.

Mount remote file system
------------------------
Note that :ref:`csiactcopy` does not include functionality to establish a conection to a remote file system.
For the moment, it is capable of copying data from one location to another. Therefore it is necessary to
mount the remote file system where the data is located on the local machine. This can be achieved using e.g.
`sshfs <https://github.com/libfuse/sshfs>`_:

.. code-block:: bash

  $ sshfs user@remote.server.com:/path/to/remote/fits/data/ /path/to/mountpoint/

The success of this operation can be verified by checking if there is a file called ``master.json``
in the mounted directory.

Copy data 
---------
After mounting the file system to a folder of your choice (``/path/to/mountpoint/`` in this case), the copying of the data can be initiated.
To monitor the progress on the screen, the hidden parameter ``debug=yes`` can be used. To increase the chattiness of the output,
the hidden parameter ``chatter`` can e.g. be set to 3.

.. code-block:: bash
  
  $ csiactcopy debug=yes chatter=3
  Location of remote master file [/path/to/mountpoint/master.json] 
  Name of FITS production to download [iact-fits]
  Destination path of FITS data [/path/to/local/fits/data/] 
  
The name of the FITS production must be available in the remote 
``master.json`` file. In case the names of the available productions are unkown
there are two options to determine an appropriate name:

*  Run :ref:`csiactdata` and set as ``datapath`` the mountpoint
*  Try any name, :ref:`csiactcopy` will throw a ``RuntimeError`` with the available production names if the input name is not found

Copy only a subset of the data
------------------------------
If the user is interested only in a small subset of the remote FITS production,
it is also possible to download only parts of the data. The user simply has to provide an
ASCII file containing a list of observation IDs to be downloaded (one ID per line). The hidden
parameter ``runlist`` then needs to be specified.

.. code-block:: bash
  
  $ csiactcopy runlist=myrunlist.txt debug=yes
  
Troubleshooting
---------------
In case a download was interupted (e.g. the remote file system was disconnected), the script can simply be rerun.
The index files which list the files available on the system are refreshed after the copying process.
The hidden boolean parameter ``clobber`` (default=yes) specifies if the user wishes that local content should be overwritten with
remote content. In turn, specifying ``clobber=no`` can speed up the process since remote content does not need to be copied
if it is already available on the user machine.

To verify the copying process was successful, it is recommended to run :ref:`csiactdata` after the process in order to list available productions
on the local machine. The output should now include the copied FITS production name.
