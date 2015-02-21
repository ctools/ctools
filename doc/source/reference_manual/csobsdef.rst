.. _csobsdef:

csobsdef
========

Generates an observation definition file.


Synopsis
--------

Generates an observation definition file.


General parameters
------------------

``inpnt [file]``
    Pointing definition ASCII file in colon separated values format.

``outobs [file]``
    Output Observation Definition XML file.

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.

``emin [real]``
    Lower energy limit of events (in TeV).
 	 	 
``emax [real]``
    Upper energy limit of events (in TeV).
 	 	 
``duration [real]``
    Duration of observation (in seconds).

``rad [real]``
    ROI radius (in degrees).

``(deadc = 0.95) [real]``
    Average deadtime correction factor.
 	 	 

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
    Specifies whether an existing output counts cube should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csobsdef.log) [string]``
    Log filename.


Related tools
-------------

None
