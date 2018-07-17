.. _cspull:

cspull
======

Generates pull distribution for all free model parameters.


Synopsis
--------

This script generates pull distributions for all free model parameters.
The pull is defined as the fitted model parameter value minus the true
value, divided by the parameter error. :ref:`cspull` will perform ``ntrials``
statistically independent Monte Carlo simulations followed by maximum
likelihood model fitting to derive the pull distribution for each free
model parameter. If the model fit is unbiased and the parameters are 
correct, the pull distribution should follow a Gaussian centred on 0
and with a sigma parameter of 1.

The script can generate pull distributions for unbinned, binned, stacked or
On/Off analyses. The script always starts from simulating an event list. If
the ``enumbins`` parameter is zero and ``onsrc=NONE`` an unbinned analysis will
be performed. If ``enumbins`` is positive and ``onsrc=NONE`` the script will
perform a binned analysis. To perform a stacked analysis, an observation definition
XML file with response cube information has to be specified for the ``inobs``
parameter. Finally, if the ``onsrc`` parameter is set to the name of a source in
the input model (specified via the ``inmodel`` parameter), the script will perform
an On/Off analysis using the ``REFLECTED`` method and a On region radius that is
given by the ``onrad`` parameter.

:ref:`cspull` supports multiprocessing. By default each simulation/analysis will
be performed in parallel over as many processes as the number of CPUs available on your
machine. The maximum number of parallel processes can be set by the user through the
``nthreads`` hidden parameter.

:ref:`cspull` will generate an ASCII file in comma-separated value (CSV) format,
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

``(inobs = NONE) [file]``
    Event list, counts cube, or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``onsrc [string]``
    Name of On source (only for On/Off analysis; specify ``NONE`` for other analyses).

``onrad [real]``
    Radius of On region (deg).

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation?

``(deadc = 0.98) [real]``
    Average deadtime correction factor.

``outfile [file]``
    ASCII file containing the individual pull values.

``ntrials [integer]``
    Number of samples for generating the pull distribution.

``ra [real]``
    Right Ascension of CTA pointing (J2000, in degrees).

``dec [real]``
    Declination of CTA pointing (J2000, in degrees).

``(rad = 5.0) [real]``
    ROI radius (in degrees).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``emin [real]``
    Lower energy limit of events (in TeV).

``emax [real]``
    Upper energy limit of events (in TeV).

``enumbins [integer]``
    Number of energy bins (0=unbinned).

``tmin [time]``
    Start time (UTC string, JD, MJD or MET in seconds).

``tmax [time]``
    Stop time (UTC string, JD, MJD or MET in seconds).

``(mjdref = 51544.5) [real]``
    Reference Modified Julian Day (MJD) for MET.

``npix [integer]``
    Number of pixels for binned analysis.

``binsz [real]``
    Pixel size for binned analysis (degrees/pixel).

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.

``(profile = no) [boolean]``
    Use likelihood profile method for errors?

``(seed = 1) [integer]``
    Initial random number generator seed.


Standard parameters
-------------------

``(nthreads = 0) [integer]``
    Number of parallel processes (0=use all available CPUs).

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = cspull.log) [string]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
:doc:`cstsdist`
