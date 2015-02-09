.. _ctulimit:

ctulimit
===========

Computes upper limit.


Synopsis
--------

Calculates the upper limit flux of a given source and stores the value in ascii file.
The input model has to be optimised for ctulimit to converge.

General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.
 	 	 
``inmodel [file]``
    Model XML file containing the source and background definitions.
 	 	 
``srcname [string]``
    Name of model component for which upper limit should be computed.
 	 	 
``expcube = NONE [file]``
    Exposure cube file (only needed for stacked analysis).

``psfcube = NONE [file]``
    PSF cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.
 	 	 
``outfile [file]``
    Output ASCII file name.

``(cl = 0.95) [real]``
    Confidence Level of upper limit
    
``(sigma_min = 0.0) [real]``
    Minimum boundary to start searching for upper limit value. Number of standard deviations above best fit value
    
``(sigma_max = 10.0) [real]``
    Maximum boundary to start searching for upper limit value. Number of standard deviations above best fit value  
 	 	 
``(emin = 1.0) [real])``
    Minimum energy of flux upper limit (in TeV).
 	 	 
``(emax = 100) [real]``
    Maximum energy of flux upper limit (in TeV).
 	 	 
``(tol = 1e-5) [real]``
    Computation tolerance for minimum
   
``(max_iter = 50) [integer]``
    Maximum number of iterations before throwing an exception


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

``(logfile = ctskymap.log) [string]``
    Name of log file.


Related tools
-------------

None
