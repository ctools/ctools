.. _cscaldb:

cscaldb
=======

Dumps all available instrument response functions into log file.


Synopsis
--------

This scripts dumps all instrument response functions that are available
in the calibration database into the log file. If you want the information
on the screen, add the ``debug=yes`` parameter.


General parameters
------------------

None
 	 	 

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
    Specifies whether an existing output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").


Related tools of scripts
------------------------

None
