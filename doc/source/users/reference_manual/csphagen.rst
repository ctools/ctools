.. _csphagen:

csphagen
=========

Creates files necessary to perform a classical region-based spectral analysis on IACT data, i.e., source and background PHA files, ARF and RMF files.


Synopsis
--------

This script can be used to derive from an observation or a set of observations the count spectra in a source region and in background regions, as well as the detector response (effective area, energy redistribution matrix), that can be used to perform a classical 1D spectral analysis. The output files are saved in the OGIP format normally used in X-ray astronomy (PHA, ARF, RMF). `See here <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html>`__.

The outputs are: 1) the PHA, ARF, RMF files, either separately for each observation, or stacked for a set of observations, 2) DS9 regions files listing the source and background regions for each observation; 3) a new observation definition XML file. 


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file

``exclusion [file]``
    Optional FITS file containing a WCS map in the first edu defining sky regions not to be used for background estimation (where map value > 0) 

``outroot [string]``
    Root of the file name for output PHA, ARF, RMF, XML, and ds9 reg files

``emin [real]``
    Lower energy limit (TeV) if LIN or LOG binning algorithms are used

``emax [real]``
    Upper energy limit (TeV) if LIN or LOG binning algorithms are used

``enumbins [integer]``,  i, a, 120,,, "Number of energy bins"
    Number of energy bins. At least 30 bins per decade are recommended for proper evaluation of the instrument response

``ebinalg [string]’’
    Energy binning algorithm: LOG - logarithmically spaced energy bins; LIN - linearly spaced energy bins; FILE - energy bounds retrieved from file

``ebinfile [file]``
    Name of the file containing the energy bin definition if FILE algorithm is used

``bkgmethod [string]’’
    Method for background estimation. REFLECTED: background evaluated in regions with the same shape as source region reflected w.r.t. pointing direction 

``srcshape [string]``
    Shape of the source region. CIRCLE: circular region around given position

``coordsys [string]``
    Coordinate system (CEL - celestial, GAL - galactic)

``ra [real]``
    Right Ascension of source region centre (deg)

``dec [real]``
    Declination of source region centre (deg)

``glon [real]``
    Galactic longitude of source region centre (deg)

``glat [real]``
    Galactic latitude of source region centre (deg)

``rad [real]``
    Radius of source region circle (deg

``bkgregmin [integer]``
    Minimum number of background region for REFLECTED method. If not available, observation is skipped. (Default: 2)

``maxoffset [real]``
    Maximum offset of source from camera center to accept observation

``stack [boolean]``
    Yes: Stack multiple observations into single PHA. No: produce run-by-run PHA files


Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged
     
     ``chatter = 1``: only errors will be logged
     
     ``chatter = 2``: errors and actions will be logged
     
     ``chatter = 3``: report about the task execution
     
     ``chatter = 4``: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output runlist should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csfindobs.log) [filename]``
    Log filename.



