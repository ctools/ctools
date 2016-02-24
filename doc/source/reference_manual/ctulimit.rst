.. _ctulimit:

ctulimit
========

Computes upper flux limit for a source model.


Synopsis
--------

This tool computes the upper flux limit for a specific source model. Except
of the node function, all spectral models are supported. Starting from the
maximum likelihood model parameters, the tool finds the model flux that leads
to a decrease of the likelihood that corresponds to a given confidence level.
By default a confidence level of 95% is used, but this level can be adjusted
using the hidden ``confidence`` parameter.

ctulimit writes the differential upper flux limit at a given reference 
energy (specified by the hidden parameter ``eref``) and the integrated 
upper flux limit (specified by the hidden parameters ``emin`` and ``emax``)
into the log file.



General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.
 	 	 
``inmodel [file]``
    Input model XML file.
 	 	 
``srcname [string]``
    Name of source model for which the upper flux limit should be computed.
 	 	 
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
 	 	 
``(eref = 1.0) [real])``
    Reference energy for differential limit (in TeV).
 	 	 
``(emin = 1.0) [real])``
    Minimum energy for integral flux limit (in TeV).
 	 	 
``(emax = 100) [real]``
    Maximum energy for integral flux limit (in TeV).
 	 	 
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
    Specifies whether an existing output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctulimit.log) [string]``
    Name of log file.


Related tools
-------------

:ref:`ctlike`
:ref:`cterror`
:ref:`ctbutterfly`
