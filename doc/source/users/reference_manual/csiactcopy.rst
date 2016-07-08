.. _csiactcopy:

csiactcopy
===========

Copies IACT data from one location to another


Synopsis
--------

This script allows to download IACT data from a remote machine. The remote file system
has to be mounted before running the script. For mounting, e.g. sshfs can be used:

.. code-block:: bash

  $ sshfs user@remote.server.com:/path/to/remote/fits /path/to/mountpoint/

The script can optionally take an ASCII file with a list of observation IDs as input.
In this way the user can copy only a specific subset of the data. Index files and data
structure is updated accordingly by the script. 

In case the copying procedure failed due to e.g. a broken connection, the script can be executed
again using the parameter ``clobber=no``. Thus files that were already copied don't get
recopied and overwritten which saves time.

In order to monitor the progress on the screen, the script can be executed with the
hidden parameter ``debug=yes``


General parameters
------------------

``remote_master [file]``
    Location of remote master file

``prodname [string]``
    Name of FITS production to download
    
``outpath [file]``
    Destination path of FITS data

``(runlist = NONE) [file]``
    List of observation IDs

    
Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether existing output files should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csiactdload.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`csiactdata`