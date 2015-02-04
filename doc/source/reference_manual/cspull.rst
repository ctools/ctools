cspull
======

Generates pull distribution for all model parameters.


Synopsis
--------

Generates pull distributions for all parameters in a model.
The pull is defined by the fitted value minus the simulated true value,
devided by the statistical uncertainty.
For unbiased parameter estimates, the pull distribution should follow a
Gaussian centred on 0 and a sigma parameter of 1.


General parameters
------------------

``inmodel [file]``
    XML file that describes the astrophysical sources and the instrumental
    background.
 	 	 
``outfile [file]``
    ASCII file containing the individual pull values.
 	 	 
``ntrials [integer]``
    Number of samples for generating the pull distribution.
 	 	 
``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation.

``ra [double]``
    Right Ascension of CTA pointing (J2000, in degrees).
 	 	 
``dec [double]``
    Declination of CTA pointing (J2000, in degrees).
 	 	 
``emin [double]``
    Lower energy limit of events (in TeV).
 	 	 
``emax [double]``
    Upper energy limit of events (in TeV).
 	 	 
``enumbins [integer]``
    Number of energy bins (0=unbinned).
 	 	 
``duration [double]``
    Duration of observation (in seconds).
 	 	 
``(deadc = 0.95) [double]``
    Average deadtime correction factor.
 	 	 
``(rad = 5.0) [double]``
    ROI radius (in degrees).

``(pattern = single) [string]``
    Pattern for pointing simulation (single/four).

``(offset = 1.5) [double]``
    Observation pattern offset (in degrees).
 	 	 
``(npix = 200) [integer]``
    Number of pixels for binned analysis.
 	 	 
``(binsz = 0.05) [double]``
    Pixel size for binned analysis.


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


Related tools
-------------

None
