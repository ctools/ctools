.. _cterror:

cterror
=======

Computes parameter errors for a specific sky model component using
a likelihood profile method.


Synopsis
--------

This tool computes the parameter errors for a specific sky model using
a likelihood profile method.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.
 	 	 
``inmodel [file]``
    Model XML file containing the source and background definitions.
 	 	 
``outmodel [file]``
    Model XML file with updated error information.
 	 	 
``srcname [string]``
    Name of model component for which upper limit should be computed.
 	 	 
``expcube [file]``
    Exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Background cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.
 	 	 
``(confidence = 0.68) [real]``
    Confidence level for error computation.
    
``(tol = 1e-3) [real]``
    Computation tolerance.
   
``(max_iter = 50) [integer]``
    Maximum number of iterations before stopping the likelihood
    profil computations.


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

``(logfile = cterror.log) [string]``
    Name of log file.


Related tools
-------------

<ctulimit>
