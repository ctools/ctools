.. _ctmodel:

ctmodel
=======

Computes model counts cube.


Synopsis
--------

This tool computes a model counts cube based on an XML model.
The counts cube is a 3-dimensional data cube spanned by Right Ascension or
Galactic longitude, Declination or Galactic latitude, and energy.


General parameters
------------------

``incube = NONE [file]``
    Input counts cube. If a filename is given the counts cube will be
    used to define the dimensions of the model cube. The content of
    the counts cube is not relevant.

``inmodel = $CTOOLS/share/models/crab.xml [string]``
    Source model XML file.

``outcube = modcube.fits [file]``
    Output model counts cube.
 	 	 
``(edisp = no) [boolean]``
    Apply energy dispersion to response computation.

``inobs = NONE [file]``
    Input event list, counts cube or observation definition XML file.
    This information is needed to determine the observation information
    (number of observations, duration, liftime, pointing, eventually
    observation specific instrument response functions, etc.).

``expcube = NONE [file]``
    Exposure cube file (only needed for stacked analysis).

``psfcube = NONE [file]``
    PSF cube file (only needed for stacked analysis).

``caldb = dummy [string]``
    Calibration database.
 	 	 
``irf = cta_dummy_irf [string]``
    Instrument response function.
 	 	 
``ra = 83.63 [real]``
    Right Ascension of pointing (J2000, in degrees).
 	 	 
``dec = 22.01 [real]``
    Declination of pointing (J2000, in degrees).

``rad = 5.0 [real]``
    Radius of field of view (in degrees).
 	 	 
``tmin = 0.0 [real]``
    CTA mission elapsed start time (in seconds).
 	 	 
``tmax = 1800.0 [real]``
    CTA mission elapsed stop time (in seconds).
 	 	 
``(deadc = 0.95) [real]``
    Average deadtime correction factor.

``(ebinalg = LOG) <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``emin = 0.1 [real]``
    Lower energy value for first energy bin (in TeV).
 	 	 
``emax = 100 [real]``
    Upper energy value for last energy bin (in TeV).
 	 	 
``enumbins = 20 [integer]``
    Number of energy bins.
 	 	 
``(ebinfile = NONE) [file]``
    Name of the file containing the energy bin definition.

``(usepnt = no) [boolean]``
    Use pointing instead of xref/yref parameters?
 	 	 
``nxpix = 200 [integer]``
    Size of the Right Ascension / Galactic longitude axis (in pixels).
 	 	 
``nypix = 200 [integer]``
    Size of the Declination / Galactic latitude axis (in pixels).
 	 	 
``binsz = 0.02 [real]``
    Pixel size (in degrees/pixel).
 	 	 
``coordsys = CEL <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``xref = 83.63 [real]``
    Right Ascension / Galactic longitude of image centre (J2000, in degrees).
 	 	 
``yref = 22.01 [real]``
    Declination / Galactic latitude of image centre (J2000, in degrees).
 	 	 
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

``(logfile = ctmodel.log) [string]``
    Name of log file.


Related tools
-------------

None
