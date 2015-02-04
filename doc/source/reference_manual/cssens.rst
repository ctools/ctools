cssens
======

Computes the CTA sensitivity for a number of energy bins.


Synopsis
--------

Computes the CTA sensitivity for a number of energy bins using :doc:`ctlike`.
Crab spectra are fitted in narrow energy bins to simulated data, and the
flux level is determined that leads to a particular significance.


General parameters
------------------

``outfile [file]``
    ASCII file containing the sensitivity values.
 	 	 
``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.
 	 	 
``type <point|gauss|shell|disk> [string]``
    Source model type.
 	 	 
``offset [double]``
    Source offset angle (in degrees).
 	 	 
``bkg [string]``
    Background model file function (none=power law for subarray E).
 	 	 
``duration [double]``
    Effective exposure time (in seconds).
 	 	 
``rad [double]``
    Radius of ROI (in degrees).
 	 	 
``(enumbins = 0) [integer]``
    Number of energy bins (0=unbinned).
 	 	 
``(npix = 200) [integer]``
    Number of pixels for binned analysis.
 	 	 
``(binsz = 0.05) [double]``
    Pixel size for binned analysis.
 	 	 
``(sigma = 5.0) [double]``
    Significance threshold.
 	 	 
``(ts_use = yes) [boolean]``
    Use TS for signifiance estimation?
 	 	 
``(index = -2.48) [double]``
    Spectral index for source model.
 	 	 
``(radius = 0.1) [double]``
    Extended source model radius (in degrees).
 	 	 
``(width = 0.05) [double]``
    Extended source model width (in degrees).
 	 	 
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


Related tools
-------------

None
