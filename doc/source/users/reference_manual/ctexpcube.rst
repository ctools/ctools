.. _ctexpcube:

ctexpcube
=========

Generate exposure cube for a stacked analysis.


Synopsis
--------

This tool generates an exposure cube for a stacked analysis. An exposure
cube is a 3-dimensional cube spanned by Right Ascension or Galactic longitude,
Declination or Galactic latitude, and energy, that gives the exposure 
as function of true sky direction and energy. The energy binning of the cube 
may be either linear, logarithmic, or custom defined using an input file.

:ref:`ctexpcube` requires on input the event list or observation definition file
that has been used in the generation of the counts cube using :doc:`ctbin`.

The binning of an exposure cube does not need to correspond to the binning
of a counts cube, as exposure values will in any case be determined by linear
interpolation. For the user's convenience, however, the counts cube can be
specified as ``incube`` parameter which will instruct ctexpcube to extract
the exposure cube definition (such as sky coordinates and projection, number
of pixels, pixel scale, energy binning) from the counts cube. Note, however,
that the exposure cube is defined in true sky coordinates and energy while
the counts cube is defined in measured sky coordinates and energy. Consequently,
the sky area and energy range covered by the exposure cube should be slightly
larger than that of the counts cube to accommodate for spill over of events
due to the point spread function and energy dispersion. By computing the
exposure cube on the same grid as the counts cube, the spill over of events
from sources at the edge of cube will not be handled correctly.

Note that there should at least be 25 bins per energy decade in the exposure
cube to assure that the energy dependence of the exposure is sufficiently
well sampled (in particular at low energies).

:ref:`ctexpcube` generates an exposure cube FITS file comprising three extensions.
The primary extension contains a 3-dimensional image that contains the 
exposure values. The values are in units of cm**2 s. The next extension named
``ENERGIES`` contains a binary table that defines the energies of the exposure
cube. The last extension named ``GTI`` contains a binary table that provides the
Good Time Intervals of all observations that have been used for the computation
of the exposure cube.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``incube [file]``
    Input counts cube file to extract exposure cube definition.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``outcube [file]``
    Output exposure cube file.

``ebinalg <FILE|LIN|LOG|POW> [string]``
    Algorithm for defining energy bins. For ``FILE``, the energy bins are defined
    in a FITS file that is specified by the ``ebinfile`` parameter, for ``LIN``
    ``LOG`` and ``POW`` there will be ``enumbins`` energy bins spaced linearly,
    logarithmically, or following a power law between ``emin`` and ``emax``,
    respectively. For ``POW``, the parameter ``ebingamma`` specifies the slope
    of the power law.

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

``ebingamma [real]``
    Exponent of the power law for ``POW`` energy binning. An exponent of 1.0
    corresponds to a logarithmic energy binning.

``(addbounds = no) [boolean]``
    Add energies to the exposure cube at the observation energy boundaries?

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

``(publish = no) [boolean]``
    Specifies whether the exposure cube should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing output exposure cube file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = ctexpcube.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctbin`
