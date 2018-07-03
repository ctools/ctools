.. _cslightcrv:

cslightcrv
==========

Computes a light curve for a given source.


Synopsis
--------

This script computes a light curve by performing a maximum likelihood fit
using :doc:`ctlike` in a series of time bins. The time bins can be either
specified in an ASCII file, as an interval divided into equally sized time
bins, or can be taken from the Good Time Intervals of the observation(s).
The format of the ASCII file is one row per time bin, each specifying the
start of stop value of the bin, separated by a whitespace. The times are
given in Modified Julian Days (MJD).

:ref:`cslightcrv` can perform the light curve fitting in unbinned, stacked or
On/Off analysis mode; on input the script always requires unbinned data, i.e.
event list(s), since only event data retain time information for light curve
fitting. To select unbinned analysis, specify ``method=3D`` and ``enumbins=0``;
for stacked analysis specify ``method=3D`` and any positive number for ``enumbins``.
To perform the analysis in On/Off model, specify ``method=ONOFF``.

For a stacked analysis, :ref:`cslightcrv` will query the parameters that define
the counts cube (``coordsys``, ``proj``, ``xref``, ``yref``, ``nxpix``, ``nypix``,
``binsz``) while for an On/Off analysis the script will query the parameters that
are necessary to run internally the :ref:`csphagen` script (``inexclusion``, ``coordsys``,
``xref``, ``yref``, ``srcshape``, ``rad``, ``bkgmethod``, ``bkgregmin``, ``maxoffset``,
``etruemin``, ``etruemax``, ``etruebins``).

:ref:`cslightcrv` writes the fitted model parameters and their statistical errors
in a FITS file. In addition, it computes for each time bin the statistical 
significance of the detection, expressed by the Test Statistics, and the 
upper flux limit.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for sensitivity computation.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

``(inexclusion = NONE) [file]``
    Optional FITS file containing a WCS map in the first hdu that defines sky
    regions not to be used for background estimation in On/Off analysis (where
    map value != 0).

``(edisp = no) [boolean]``
    Applies energy dispersion to response computation (for ``3D`` analysis only,
    energy dispersion is always taken into account in On/Off analysis).

``outfile [file]``
    Name of the light curve output file.

``tbinalg <FILE|LIN|GTI> [string]``
    Algorithm for defining time bins.

``tmin [time]``
    Lightcurve start time (UTC string, JD, MJD or MET in seconds).

``tmax [time]``
    Lightcurve stop time (UTC string, JD, MJD or MET in seconds).

``(mjdref = 51544.5) [real]``
    Reference Modified Julian Day (MJD) for Mission Elapsed Time (MET).

``tbins [integer]``
    Number of time bins.

``tbinfile [file]``
    File defining the time binning.

``method  <3D|ONOFF> [string]``
    Selects between ``3D`` analysis (3D spatial/energy likelihood) and ``ONOFF``
    analysis (1D likelihood with background from Off regions).

``emin [real]``
    Lower energy limit of events (in TeV).

``emax [real]``
    Upper energy limit of events (in TeV).

``enumbins [integer]``
    Number of energy bins per light curve bin (0=unbinned for ``3D`` analysis only).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of image centre for ``3D`` analysis or
    source region centre for On/Off analysis (J2000, in degrees).

``yref [real]``
    Declination / Galactic latitude of image centre for ``3D`` analysis or
    source region centre for On/Off analysis (J2000, in degrees).

``nxpix [integer]``
    Size of the Right Ascension / Galactic longitude axis for ``3D`` analysis (in pixels).

``nypix [integer]``
    Size of the Declination / Galactic latitude axis for ``3D`` analysis (in pixels).

``binsz [real]``
    Pixel size for ``3D`` analysis (in degrees/pixel).

``(srcshape = CIRCLE) [string]``
    Shape of the source region for On/Off analysis.
    ``CIRCLE``: circular region around given position.

``rad [real]``
    Radius of source region circle for On/Off analysis (deg)

``(bkgmethod = REFLECTED) [string]``
    Method for background estimation in On/Off analysis.
    ``REFLECTED:`` background evaluated in regions with the same shape as
    source region reflected w.r.t. pointing direction for each observation.

``(bkgregmin = 2) [integer]``
    Minimum number of background regions that are required for an observation in
    ``ONOFF`` analysis. If this number of background regions is not available the
    observation is skipped.

``(maxoffset = 4.0) [real]``
    Maximum offset in degrees of source from camera center to accept the
    observation for On/Off analysis.

``(etruemin = 0.01) [real]``
    Minimum true energy to evaluate instrumental response in On/Off analysis (TeV).

``(etruemax = 0.01) [real]``
    Maximum true energy to evaluate instrumental response in On/Off analysis (TeV).

``(etruebins = 30) [integer]``
    Number of bins per decade for true energy bins to evaluate instrumental
    response in On/Off analysis.

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.

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

``(nthreads = 0) [integer]``
    Number of parallel processes (0=use all available CPUs).

``(publish = no) [boolean]``
    Specifies whether the light curve should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing light curve output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = cslightcrv.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
:doc:`csphagen`
