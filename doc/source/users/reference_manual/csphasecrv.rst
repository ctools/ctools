.. _csphasecrv:

csphasecrv
==========

Computes phase dependent spectra for a given source.


Synopsis
--------

This script computes spectra by performing a maximum likelihood fit using
:doc:`ctlike` in a series of phase bins for periodic sources. The input event
list must have a ``PHASE`` column containing the phase values, that can be
assigned using :doc:`ctphase`. The phase bins can be either specified in an ASCII
file or as an interval divided into equally sized phase bins. The format of the
ASCII file is one row per phase bin, each specifying the start of stop value of
the phase bin, separated by a whitespace. The phase goes from 0.0 to 1.0.

:ref:`csphasecrv` can perform the phase curve fitting in unbinned, stacked or
On/Off analysis mode; on input the script always requires unbinned data, i.e.
event list(s), since only event data retain phase information for phase curve
fitting. To select unbinned analysis, specify ``method=3D`` and ``enumbins=0``;
for stacked analysis specify ``method=3D`` and any positive number for ``enumbins``.
To perform the analysis in On/Off model, specify ``method=ONOFF``.

For a stacked analysis, :ref:`csphasecrv` will query the parameters that define
the counts cube (``coordsys``, ``proj``, ``xref``, ``yref``, ``nxpix``, ``nypix``,
``binsz``) while for an On/Off analysis the script will query the parameters that
are necessary to run internally the :ref:`csphagen` script (``inexclusion``, ``coordsys``,
``xref``, ``yref``, ``srcshape``, ``rad``, ``bkgmethod``, ``bkgregmin``, ``maxoffset``,
``etruemin``, ``etruemax``, ``etruebins``).

:ref:`csphasecrv` supports multiprocessing. By default the analysis in each phase bin
will be performed in parallel over as many processes as the number of CPUs
available on your machine. The maximum number of parallel processes can be set
by the user through the ``nthreads`` hidden parameter.

The script writes the fit results into a FITS file. The script also produces one
XML file per phase bin that contains the best-fit model for that bin.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``srcname [string]``
    Name of the periodic source in the source model XML file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

``(inexclusion = NONE) [file]``
    Optional FITS file containing a WCS map that defines sky regions
    not to be used for background estimation (where map value !=
    0). If the file contains multiple extensions the user may specify
    which one to use. Otherwise, the extention ``EXCLUSION`` will be
    used if available, or else the primary extension will be used.

``(edisp = no) [boolean]``
    Applies energy dispersion to response computation (for ``3D`` analysis only,
    energy dispersion is always taken into account in ``ONOFF`` analysis).

``outfile [file]``
    Name of the XML output file. The phase interval will be automatically
    appended to the name.

``phbinalg <FILE|LIN> [string]``
    Algorithm for defining phase bins.

``phbins [integer]``
    Number of phase bins.

``phbinfile [file]``
    File defining the phase binning.

``method <3D|ONOFF> [string]``
    Selects between 3D analysis (3D spatial/energy likelihood) and ONOFF
    analysis (1D likelihood with background from Off regions).

``emin [real]``
    Lower energy limit of events (in TeV).

``emax [real]``
    Upper energy limit of events (in TeV).

``enumbins [integer]``
    Number of energy bins per phase bin (0=unbinned for ``3D`` analysis only).

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
    On/Off analysis. If this number of background regions is not available the
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


Standard parameters
-------------------

``(nthreads = 0) [integer]``
    Number of parallel processes (0=use all available CPUs).

``(publish = no) [boolean]``
    Specifies whether the phase curve should be published on VO Hub.

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

``(logfile = csphasecrv.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctphase`
:doc:`ctlike`
