.. _csspec:

csspec
======

Computes spectrum for a given source.


Synopsis
--------

Computes spectral points for a given source.


General parameters
------------------

``(inobs = NONE) [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input source model XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for sensitivity computation.

``outfile [file]``
    Name of the source spectrum output file.

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

``emin [real]``
    Lower energy limit of events (in TeV).
 	 	 
``emax [real]``
    Upper energy limit of events (in TeV).
 	 	 
``enumbins [integer]``
    Number of energy bins (0=unbinned).
 	 	 
``(ebinalg = LOG) <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``binned [boolean]``
    Specifies whether a binned computation should be used.

``(nebins = 5) [integer]``
    Number of bins per spectral point for binned analysis.

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``proj <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of image centre (J2000, in degrees).
 	 	 
``yref [real]``
    Declination / Galactic latitude of image centre (J2000, in degrees).

``nxpix [integer]``
    Size of the Right Ascension / Galactic longitude axis (in pixels).
 	 	 
``nypix [integer]``
    Size of the Declination / Galactic latitude axis (in pixels).
 	 	 
``binsz [real]``
    Pixel size (in degrees/pixel).
 	 	 
``(calc_ts = yes) [boolean]``
    Compute TS for each spectral point.

``(calc_ulim = yes) [boolean]``
    Compute upper limit for each spectral point.

``(fix_srcs = yes) [boolean]``
    Fix other sky model parameters.

``(fix_bkg = no) [boolean]``
    Fix background parameters.


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

``(logfile = csspec.log) [filename]``
    Log filename.


Related tools
-------------

None
