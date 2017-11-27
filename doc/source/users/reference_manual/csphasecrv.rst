.. _csphasecrv:

csphasecrv
==========

Computes phase dependent spectra for a given source.


Synopsis
--------

This script computes spectra by performing a maximum likelihood fit using
:doc:`ctlike` in a series of phase bins for periodic sources. The input event
list must have a ``PHASE`` column containing the phase values, that can be
assigned using :doc:`ctphase`. The phase bins can be
either specified in an ASCII file or as an interval divided into equally sized
phase bins. The format of the ASCII file is one row per phase bin, each
specifying the start of stop value of the phase bin, separated by a whitespace.
The phase goes from 0.0 to 1.0.

The script writes the fit results into a FITS file. The script
also produces one XML file per phase bin that contains the best-fit model for
that bin.


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
    Optional FITS file containing a WCS map in the first hdu that defines sky
    regions not to be used for background estimation in On/Off analysis (where
    map value != 0).

``(edisp = no) [boolean]``
    Applies energy dispersion to response computation (for Cube analysis only,
    energy dispersion is always taken into account in On/Off analysis).

``outfile [file]``
    Name of the XML output file. The phase interval will be automatically
    appended to the name.

``phbinalg <FILE|LIN> [string]``
    Algorithm for defining phase bins.

``phbins [integer]``
    Number of phase bins.

``phbinfile [file]``
    File defining the phase binning.

``method [string]``
    Selects between CUBE analysis (3D spatial/energy likelihood) and ONOFF
    analysis (1D likelihood with background from Off regions).

``emin [real]``
    Lower energy limit of events (in TeV).

``emax [real]``
    Upper energy limit of events (in TeV).

``enumbins [integer]``
    Number of energy bins per phase bin (0=unbinned).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
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

``rad [real]``
    Radius of source region circle for On/Off analysis (deg)

``(bkgmethod = REFLECTED) [string]``
    Method for background estimation in On/Off analysis.
    ``REFLECTED:`` background evaluated in regions with the same shape as
    source region reflected w.r.t. pointing direction for each observation.

``(bkgregmin = 2) [integer]``
    Minimum number of background regions that are required for an observation in
    On/Off analysis. If this number of background regions is not available the observation is
    skipped.

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
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csphasecrv.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctphase`
:doc:`ctlike`
