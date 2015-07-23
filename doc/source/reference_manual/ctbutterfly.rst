.. _ctbutterfly:

ctbutterfly
===========

Computes butterfly.


Synopsis
--------

Calculates the confidence band of a spectral model, taking into account the
covariance matrix of a likelihood fit.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.
 	 	 
``inmodel [file]``
    Model XML file containing the source and background definitions.
 	 	 
``srcname [string]``
    Name of model component for which butterfly should be computed.
 	 	 
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

``(edisp = no) [boolean]``
    Applies energy dispersion to response computation.
 	 	 
``outfile [file]``
    Output butterfly ASCII file name.
 	 	 
``(matrix = "NONE") [file]``
    Input covariance matrix file (not used)

``(ebinalg = "LOG") <FILE|LIN|LOG> [string]``
    Butterfly energy binning algorithm.
 	 	 
``emin [real]``
    Minimum energy of butterfly (in TeV).
 	 	 
``emax [real]``
    Maximum energy of butterfly (in TeV).
 	 	 
``(enumbins = 100) [string]``
    Number of energy bins of butterfly.
 	 	 
``(ebinfile = "NONE") [string]``
    File for energy binning algorithm.


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
