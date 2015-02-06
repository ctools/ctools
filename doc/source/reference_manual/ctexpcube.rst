.. _ctexpcube:

ctexpcube
=========

Generate an exposure cube for binned maximum likelihood analysis.


Synopsis
--------

This tool generates an exposure cube for use in a binned maximum likelihood
analysis.


General parameters
------------------

``inobs = NONE [file]``
    Input event list or observation definition XML file.

``incube = NONE [file]``
    Counts cube for PSF cube definition.

``outcube = expcube.fits [file]``
    Output exposure cube file.

``caldb = dummy [string]``
    Calibration database.

``irf = cta_dummy_irf [string]``
    Response function.

``(edisp = no) [boolean]``
    Apply energy dispersion for response computation?

``ebinalg = LOG <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``emin = 0.1 [real]``
    Lower energy value for first energy bin (in TeV).
 	 	 
``emax = 100.0 [real]``
    Upper energy value for last energy bin (in TeV).
 	 	 
``enumbins =20 [integer]``
    Number of energy bins.
 	 	 
``ebinfile = NONE [file]``
    Name of the file containing the energy bin definition.
 	 	 
``(usepnt = no) [boolean]``
    Use CTA pointing direction for cube centre instead of xref/yref parameters?
 	 	 
``nxpix = 200 [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.
 	 	 
``nypix = 200 [integer]``
    Number of cube bins in Declination or Galactic latitude.
 	 	 
``binsz = 0.02 [real]``
    Cube bin size (in degrees/pixel).
 	 	 
``coordsys = CEL <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``xref = 83.63 [real]``
    Right Ascension / Galactic longitude of cube centre (J2000, in degrees).
 	 	 
``yref = 22.01 [real]``
    Declination / Galactic latitude of cube centre (J2000, in degrees).
 	 	 
``(axisrot = 0.0) [real]``
    Rotation angle of image axes (in degrees).
 	 	 
``proj = CAR <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.
 	 	 

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

``(logfile = ctexpcube.log) [string]``
    Name of log file.


Related tools
-------------

None
