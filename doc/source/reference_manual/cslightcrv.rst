.. _cslightcrv:

cslightcrv
==========

Computes a lightcurve for a given source.


Synopsis
--------

Computes a lightcurve by performing a maximum likekihood fit in a series of 
time bins. The time bins can be either specified in an ASCII file, as an
interval divided into equally sized time bins, or can be taken from the Good 
Time Intervals of the observation(s). The format of the ASCII file is one
row per time bin, each specifying the start of stop value of the bin, separated
by a whitespace. The times are given in Modified Julian Days (MJD). 
cslightcrv writes the fitted model parameters and their statistical errors 
in a FITS file. In addition, it computes for each time bin the statistical 
significance of the detection, expressed by the Test Statistics, and the 
upper flux limit.


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

``tmin [real]``
    Lightcurve start time (in MJD).

``tmax [real]``
    Lightcurve stop time (in MJD).

``tbins [integer]``
    Number of time bins.

``(tbinalg = GTI) <FILE|LIN|GTI> [string]``
    Algorithm for defining time bins.

``tbinfile [file]``
    File defining the time binning.

``binned [boolean]``
    Use binned analysis in each time bin.

``enumbins [integer]``
    Number of energy bins per light curve bin (0=unbinned).

``emin [real]``
    Lower energy limit of events (in TeV).
 	 	 
``emax [real]``
    Upper energy limit of events (in TeV).
 	 	 
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
    Compute TS value for each time bin?

``(calc_ulim = yes) [boolean]``
    Compute upper limit for each time bin?

``(fix_srcs = yes) [boolean]``
    Fix other sky model parameters?

``(fix_bkg = no) [boolean]``
    Fix background model parameters?


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

``(logfile = cslightcrv.log) [filename]``
    Log filename.


Related tools
-------------

None
