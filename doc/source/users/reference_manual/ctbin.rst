.. _ctbin:

ctbin
=====

Generate counts cube from event list(s).


Synopsis
--------

This tool creates a counts cube and fills it with events. A counts cube is 
a 3-dimensional data cube spanned by Right Ascension or Galactic longitude,
Declination or Galactic latitude, and energy. The energy binning may be either
linear, logarithmic, or custom defined using an input file. The events are 
either taken from a single event list file or from the event lists that are 
specified in an observation definition file.

In case that multiple event lists are given in an observation definition file
there are two options. By default, the tool will loop over all event lists and
fill all events into a single counts cube. If ``stack=no`` is specified, the
tool will bin each event list into a separate counts cube.

:ref:`ctbin` generates counts cube FITS file(s) comprising three extensions. The
primary extension contains a 3-dimensional image that contains the counts
cube values. The next extension named ``EBOUNDS`` contains a binary table
that defines the energy boundaries of the counts cube. The last extension
named ``GTI`` contains a binary table that defines the Good Time Intervals
of all event lists that have been filled into the counts cube.

By default, the tool will write out a single, stacked, counts cube. If ``stack=no``
is specified, the tool will write all counts cubes into a separate FITS file,
prefixing the input event list file name with the string that is specified by
the hidden ``prefix`` parameter. In addition, the tool will create an observation
definition XML file that lists all produced counts cubes, and that can be used
as input for any further analysis step.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``outobs [file]``
    Output counts cube file or observation definition XML file.

``(stack = yes) [boolean]``
    Specify whether stacking of the events in a single counts cube is requested.
    By default this hidden parameter is set to yes, and the tool produces a
    single counts cube FITS file on output. If ``no`` is specified, the tool
    will produce a counts cube for each event list that is found in the input
    observation, and will also write an observation definition XML file.

``(prefix = cntcube_) [string]``
    Prefix for output counts cube in observation definition XML file.

``ebinalg <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins. For ``FILE``, the energy bins are defined
    in a FITS file that is specified by the ``ebinfile`` parameter, for ``LIN``
    and ``LOG`` there will be ``enumbins`` energy bins spaced linearly or
    logarithmically between ``emin`` and ``emax``, respectively.

``emin [real]``
    Lower energy value for first energy bin (in TeV) if ``LIN`` or ``LOG``
    energy binning algorithms are used.

``emax [real]``
    Upper energy value for last energy bin (in TeV) if ``LIN`` or ``LOG``
    energy binning algorithms are used.

``enumbins [integer]``
    Number of energy bins if ``LIN`` or ``LOG`` energy binning algorithms are
    used.

``ebinfile [file]``
    Name of the file containing the energy binning definition if ``ebinalg=FILE``.
    You may use :ref:`csebins` to generate a file with appropriate energy binning.

``(usepnt = no) [boolean]``
    Use CTA pointing direction for cube centre instead of xref/yref parameters?

``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.

``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.

``binsz [real]``
    Cube bin size (in degrees/pixel).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of cube centre (J2000, in degrees).

``yref [real]``
    Declination / Galactic latitude of cube centre (J2000, in degrees).


Standard parameters
-------------------

``(nthreads = 0) [integer]``
    Number of parallel processes (0=use all available CPUs).

``(publish = no) [boolean]``
    Specifies whether the counts cube should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing output counts cube file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = ctbin.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctexpcube`
:doc:`ctpsfcube`
:doc:`ctbkgcube`
