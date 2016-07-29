.. _cstsdist:

cstsdist
========

Generates the TS distribution for a given source.


Synopsis
--------

This script generates the Test Statistics (TS) distribution for a given 
source by repeatedly computing the TS value for ``ntrials`` simulated data 
sets. The Test Statistics is defined as twice the log-likelihood difference
that is obtained when fitting simulated data with and without a given
source model component.

This script either performs an unbinned or a binned simulation, stacked 
analysis is not yet supported.

cstsdist will create an ASCII file in comma-separated value (CSV) format,
containing one row per TS computation. The first row is a header row providing
the column names. The following rows give the TS value, the log-likelihood 
values of the fit with or without the source, the number of observed and 
fitted events, as well as the values and errors for all fitted parameters.

From the output file, TS distribution plots can be generated using for
example the ``show_ts_distribution.py`` script in the examples folder. The
script require matplotlib for plotting.


General parameters
------------------

``(inobs = NONE) [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for Test Statistics computation.

``outfile [file]``
    ASCII file containing the TS distribution values.

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
 	 	 
``(deadc = 0.95) [real]``
    Average deadtime correction factor.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation.

``ntrials [integer]``
    Number of Monte Carlo samples.

``ra [real]``
    Right Ascension of CTA pointing (J2000, in degrees).
 	 	 
``dec [real]``
    Declination of CTA pointing (J2000, in degrees).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``proj <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.

``emin [real]``
    Lower energy limit (in TeV).
 	 	 
``emax [real]``
    Upper energy limit (in TeV).
 	 	 
``enumbins [integer]``
    Number of energy bins (0=unbinned).
 	 	 
``(tmin = 0.0) [real]``
    Start time (in seconds).
 	 	 
``tmax [real]``
    Stop time (in seconds).
 	 	 
``npix [integer]``
    Number of pixels for binned analysis.
 	 	 
``binsz [real]``
    Pixel size for binned analysis.

``(rad = 5.0) [real]``
    Radius of ROI (in degrees).

``(pattern = single) [string]``
    Pattern for pointing simulation (single/four).

``(offset = 1.5) [real]``
    Observation pattern offset (in degrees).


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
    Specifies whether an existing output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = cstsdist.log) [string]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
:doc:`cspull`

