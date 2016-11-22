.. _cspull:

cspull
======

Generates pull distribution for all free model parameters.


Synopsis
--------

This script generates pull distributions for all free model parameters.
The pull is defined as the fitted model parameter value minus the true
value, divided by the parameter error. cspull will perform ``ntrials`` 
statistically independent Monte Carlo simulations followed by maximum
likelihood model fitting to derive the pull distribution for each free
model parameter. If the model fit is unbiased and the parameters are 
correct, the pull distribution should follow a Gaussian centred on 0
and with a sigma parameter of 1.

cspull will generate an ASCII file in comma-separated value (CSV) format,
containing one row per pull. The first row is a header row providing the 
column names. The following rows give the pull results, one row per pull. 
This includes for each parameter the parameter value, error and pull. The 
maximum likelihood value and the observed and the estimated number of counts 
are also given.

From the output file, pull distribution plots can be generated using for
example the ``show_pull_histogram.py`` script in the examples folder. The
script ``show_pull_evolution.py`` in the same folder shows the evolution
of the mean and standard deviation of the pull as function of the number
of trials. Both scripts require matplotlib for plotting.


General parameters
------------------

``inobs [file]``
    Event list, counts cube, or observation definition XML file.

``inmodel [file]``
    Input model XML file.
 	 	 
``outfile [file]``
    ASCII file containing the individual pull values.
 	 	 
``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.

``(deadc = 0.95) [real]``
    Average deadtime correction factor.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation?

``(profile = no) [boolean]``
    Use likelihood profile method for errors?

``ntrials [integer]``
    Number of samples for generating the pull distribution.
 	 	 
``ra [real]``
    Right Ascension of CTA pointing (J2000, in degrees).
 	 	 
``dec [real]``
    Declination of CTA pointing (J2000, in degrees).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``proj <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.
 	 	 
``emin [real]``
    Lower energy limit of events (in TeV).
 	 	 
``emax [real]``
    Upper energy limit of events (in TeV).
 	 	 
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
    ROI radius (in degrees).

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
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = cspull.log) [string]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
:doc:`cstsdist`
