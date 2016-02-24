.. _csiactdata:

csiactdata
===========

Dumps information about the available FITS productions of IACT data on the
screen.


Synopsis
--------

This script inspects the local IACT FITS data storage and prints the names of
available and valid data productions on the screen. The parameter ``datapath``
is only queried if the environment variable $VHEFITS is not set. The path should
be set to the location of the `master index file <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/super_index/index.html>`__.
The output of this script can be used as input for further tools that use IACT
data.


General parameters
------------------

``datapath [string]``
    Path were data are located.
    
``(master_index = master.json) [string]``
    Name of master index file.
    
    
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

``(logfile = csiactdata.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`csfindobs`