.. _ctlike:

ctlike
======

Performs binned and unbinned maximum likelihood analysis of CTA data.


Synopsis
--------

Determines source model parameters, such as flux, spectral index, source 
position, and source extent from a maximum likelihood analysis of CTA data.
The analysis can be done in a binned or an unbinned formulation of the
log-likelihood function. For binned analysis, a counts cube produced by
:doc:`ctbin` is required. For unbinned analysis, an event list processed by
:doc:`ctselect` is required. Based on the input file format, ctlike
automatically selects between binned and unbinned maximum likelihood analysis.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [string]``
    Source model input XML file.
 	 	 
``expcube [file]``
    Exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Background cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrument response function.
 	 	 
``outmodel [string]``
    Source model result XML file with values and uncertainties updated by
    the maximum likelihood fit.

``(stat = POISSON) [string]``
    Fitting statistics (POISSON or GAUSSIAN; only affects binned analysis).
 	 	 
``(edisp = no) [boolean]``
    Applies energy dispersion to response computation.

``(refit = no) [boolean]``
    Performs refitting of solution after initial fit.
 	 	 
 	 	 

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

``(logfile = ctlike.log) [string]``
    Name of log file.


Related tools
-------------

None
