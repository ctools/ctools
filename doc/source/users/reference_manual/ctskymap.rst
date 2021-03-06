.. _ctskymap:

ctskymap
========

Generate sky map from event list(s).


Synopsis
--------

This tool creates a sky map from either a single event list or event lists
provided in an observation definition file. The tool will loop over all event
lists that are provided and fill all events into a single sky map. Only events
within an energy interval spanned by ``emin`` and ``emax`` are considered.

Optionally, the tool will subtract a background model from the sky map. The
background subtraction method can be selected using the ``bkgsubtract``
parameter. By default, no background model is subtracted (method ``NONE``).
If ``IRF`` is selected, the background template that are shipped with the
Instrument Response Functions will be used for background subtraction.
If ``RING`` is selected, the background is estimated from a ring with the
inner radius specified by the ``inradius`` parameter and the other radius
specified by the ``outradius`` parameter. If a background template is available,
that template will be used for weighting the counts in the ring. Otherwise, a
constant background rate will be assumed.

:ref:`ctskymap` generates a FITS file comprising a sky map as primary extension.
If a background subtraction method was selected, the FITS file will contain
the additional extensions ``BACKGROUND`` and ``SIGNIFICANCE`` that contain
the background map and a significance map, respectively. For the latter, the
significance of the signal is computed for each sky map pixel.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``caldb [string]``
    Calibration database (only required for IRF background subtraction if no
    response information is provided by ``inobs``).

``irf [string]``
    Instrument response function (only required for IRF background subtraction
    if no response information is provided by ``inobs``).

``(inmap = NONE) [file]``
    Input sky map file containing a ``COUNTS`` and an ``ACCEPTANCE`` extension.
    Such extensions are produce by a previous run of the :ref:`ctskymap` tool with
    the ``IRF`` or ``RING`` background method.

``outmap [file]``
    Output sky map file.

``emin [real]``
    Minimum energy in map (in TeV).

``emax [real]``
    Maximum energy in map (in TeV).

``(usepnt = no) [boolean]``
    Use pointing instead of xref/yref parameters?

``nxpix [integer]``
    Size of the Right Ascension / Galactic longitude axis (in pixels).

``nypix [integer]``
    Size of the Declination / Galactic latitude axis (in pixels).

``binsz [real]``
    Pixel size (in degrees/pixel).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of image centre (J2000, in degrees).

``yref [real]``
    Declination / Galactic latitude of image centre (J2000, in degrees).


Background subtraction configuration parameters
-----------------------------------------------

``bkgsubtract <NONE|IRF|RING> [string]``
    Background subtraction method.

``roiradius [real]``
    Source region radius for ``RING`` subtraction (in degrees).

``inradius [real]``
    Inner background ring radius for ``RING`` subtraction (in degrees).

``outradius [real]``
    Outer background ring radius for ``RING`` subtraction (in degrees).

``iterations [integer]``
    Number of iterations for the automatic computation of exclusion regions.
    Sky map pixels with detection significance above the value specified by the
    ``threshold`` parameter will be defined as exclusion regions, and the ``RING``
    background will be recomputed based on these exclusion regions. Since this
    will likely change the significance of the pixels, the procedure can be
    iteratively repeated, and the ``iterations`` parameter specifies how often
    the procedure should be repeated.

``threshold [real]``
    Significance threshold above which pixels should be excluded for the
    background computation in the ``RING`` background method.

``(inexclusion = NONE) [file]``
    Exclusion region file as either a FITS map or DS9 region file.

``(usefft = yes) [boolean]``
    Specifies whether a FFT should be used for ring computation. If set to
    ``no`` the ring background will be computed directly. Direct computation is
    more accurate since the precise angular distances between pixels are taken
    into account, but FFT is faster, but less accurate since it assumes
    Euclidean distances between pixels.


Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the sky map should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing output sky map file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = ctskymap.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
