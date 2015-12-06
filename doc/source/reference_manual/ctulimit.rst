.. _ctulimit:

ctulimit
===========

Computes upper limit for a specific sky model component.


Synopsis
--------

Computes the upper limit flux of a given source.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.
 	 	 
``inmodel [file]``
    Input model XML file.
 	 	 
``srcname [string]``
    Name of model component for which upper limit should be computed.
 	 	 
``expcube [file]``
    Input exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    Input PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Input background cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.
 	 	 
``(edisp = no) [boolean]``
    Applies energy dispersion to response computation.
 	 	 
``(confidence = 0.95) [real]``
    Confidence level of upper limit.
    
``(sigma_min = 0.0) [real]``
    Minimum boundary to start searching for upper limit value.
    Number of standard deviations above best fit value
    
``(sigma_max = 10.0) [real]``
    Maximum boundary to start searching for upper limit value.
    Number of standard deviations above best fit value  
 	 	 
``(emin = 1.0) [real])``
    Minimum energy of flux limits (in TeV).
 	 	 
``(emax = 100) [real]``
    Maximum energy of flux limits (in TeV).
 	 	 
``(tol = 1e-5) [real]``
    Computation tolerance.
   
``(max_iter = 50) [integer]``
    Maximum number of iterations before stopping the upper
    limit computations.


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

``(logfile = ctulimit.log) [string]``
    Name of log file.


Related tools
-------------

:ref:`cterror`