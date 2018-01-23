.. _csphagen:

csphagen
========

Creates files necessary to perform a region-based spectral On/Off analysis.


Synopsis
--------

This script can be used to derive from an observation or a set of observations
the count spectra in a source region and in background regions, as well as the
detector response (effective area, energy redistribution matrix), that can be
used to perform a classical 1D spectral analysis. Regions can be placed
automatically or by hand. The output files are saved in the OGIP format normally
used in X-ray astronomy (PHA, ARF, RMF).
`See here <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html>`__.

The outputs are:

1) the PHA, ARF, RMF files, either separately for each observation, or stacked
   for a set of observations;
2) DS9 regions files listing the source and background regions for each
   observation;
3) a new observation definition XML file.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``(inexclusion = NONE) [file]``
    Optional FITS file containing a WCS map in the first hdu that defines sky
    regions not to be used for background estimation (where map value != 0).

``outobs [string]``
    Output observation definition XML file.

``(prefix = onoff) [string]``
    Prefix of the file name for output PHA, ARF, RMF, XML, and DS9 region files.

``emin [real]``
    Lower energy limit (TeV) if ``LIN`` or ``LOG`` binning algorithms are used.

``emax [real]``
    Upper energy limit (TeV) if ``LIN`` or ``LOG`` binning algorithms are used.

``enumbins [integer]``
    Number of energy bins. At least 30 bins per decade are recommended for
    proper evaluation of the instrument response.

``ebinalg <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins (``FILE``: energy bounds retrieved from
    file, see ``ebinfile`` parameter; ``LIN``: linearly spaced energy bins;
    ``LOG``: logarithmically spaced energy bins).

``ebinfile [file]``
    Name of the file containing the energy bin definition if ``FILE`` algorithm
    is used.

``(srcshape = CIRCLE) <CIRCLE> [string]``
    Shape of the source region (``CIRCLE``: circular region around given position).

``coordsys <CEL|GAL> [string]``
    Coordinate system (``CEL``: celestial; ``GAL``: galactic).

``ra [real]``
    Right Ascension of source region centre (deg).

``dec [real]``
    Declination of source region centre (deg).

``glon [real]``
    Galactic longitude of source region centre (deg).

``glat [real]``
    Galactic latitude of source region centre (deg).

``rad [real]``
    Radius of source region circle (deg).

``srcregfile [file]``
    Source region file (ds9 or FITS WCS map).

``bkgmethod <REFLECTED|CUSTOM> [string]``
    Method for background estimation:

    - ``REFLECTED``: background is evaluated in regions with the same shape as
      the source region reflected w.r.t. pointing direction for each observation

    - ``CUSTOM``: background is evaluated in regions specified by user. For an
      event list or a single observation in the observation definition XML file
      a region file will be queried (see ``bkgregfile`` parameter). For multiple
      observations specified in the observation definition XML file the name of
      the region file will be extracted from the ``OffRegions`` parameter that
      needs to be specified for each observation in the observation definition
      XML file. Off region files can be either ds9 region files or FITS WCS maps.

``bkgregfile [file]``
    Background regions file (ds9 or FITS WCS map).

``(bkgregmin = 2) [integer]``
    Minimum number of background regions that are required for an observation.
    If this number of background regions is not available the observation is
    skipped.

``(maxoffset = 4.0) [real]``
    Maximum offset in degrees of source from camera center to accept the
    observation.

``stack [boolean]``
    Specifies whether multiple observations should be stacked (``yes``) or
    whether run-wise PHA, ARF and RMF files should be produced (``no``).

``(etruemin = 0.01) [real]``
    Minimum true energy (TeV).

``(etruemax = 0.01) [real]``
    Maximum true energy (TeV).

``(etruebins = 30) [integer]``
    Number of bins per decade for true energy bins.


Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing output runlist should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = csphagen.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None

