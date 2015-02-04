.. _cstsdist:

cstsdist
========

Generates the TS distribution for a particular model.


Synopsis
--------

Generates the Test Statistics (TS) distribution for a particular model based on
Monte-Carlo simulations.


General parameters
------------------

``(outfile = ts.dat) [file]``
    ASCII file containing the TS distribution values.
 	 	 
``ntrials [integer]``
    Number of Monte Carlo samples.
 	 	 
``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.
 	 	 
``type <point|gauss|shell|disk> [string]``
    Source model type.
 	 	 
``(index = -2.48) [double]``
    Spectral index for source model.
 	 	 
``offset [double]``
    Source offset angle (in degrees).
 	 	 
``bkg [string]``
    Background model file function (none=power law for subarray E).
 	 	 
``emin [double]``
    Lower energy limit (in TeV).
 	 	 
``emax [double]``
    Upper energy limit (in TeV).
 	 	 
``enumbins [integer]``
    Number of energy bins (0=unbinned).
 	 	 
``duration [double]``
    Effective exposure time (in seconds).
 	 	 
``(rad = 5.0) [double]``
    Radius of ROI (in degrees).
 	 	 
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
