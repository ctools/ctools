.. _cssens:

cssens
======

Computes the CTA sensitivity for a number of energy bins.


Synopsis
--------

Computes the CTA sensitivity for a number of energy bins using :doc:`ctlike`.
Crab spectra are fitted in narrow energy bins to simulated data, and the
flux level is determined that leads to a particular significance.
The significance is estimated using the Test Statistic value.
The simplified assumption is made that the significance (in Gaussian
sigma) is the square root of the Test Statistic.


General parameters
------------------

``(inobs = NONE) [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input source model XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for sensitivity computation.

``(offset = 0.0) [real]``
    Offset angle of source in field of view (in degrees).	 

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.

``(deadc = 0.95) [real]``
    Average deadtime correction factor.
 	 	 
``outfile [file]``
    ASCII file containing the sensitivity values.
 	 	 
``type <point|gauss|shell|disk> [string]``
    Source model type.
 	 	 
``offset [real]``
    Source offset angle (in degrees).
 	 	 
``duration [real]``
    Effective exposure time (in seconds).
 	 	 
``rad [real]``
    Radius of Region of Interest (RoI) (in degrees).
 	 	 
``emin [real]``
    Lower energy limit for sensitivity computation (in TeV).
 	 	 
``emax [real]``
    Upper energy limit for sensitivity computation (in TeV).

``bins [integer]``
    Number of energy bins for sensitivity computation.
 	 	 
``(enumbins = 0) [integer]``
    Number of energy bins (0=unbinned).
 	 	 
``(npix = 200) [integer]``
    Number of pixels for binned analysis.
 	 	 
``(binsz = 0.05) [real]``
    Pixel size for binned analysis.
 	 	 
``(sigma = 5.0) [real]``
    Significance threshold.
 	 	 
``(max_iter = 50) [integer]``
    Maximum number of iterations.
 	 	 
``(num_avg = 3) [integer]``
    Number of iterations for sliding average.


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

``(logfile = cssens.log) [filename]``
    Log filename.


Related tools
-------------

None
