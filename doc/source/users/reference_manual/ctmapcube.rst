.. _ctmapcube:

ctmapcube
=========

Generate a map cube from an input model.


Synopsis
--------

This tool generates a map cube from all sky model components contained in
an input model. A map cube is a series of sky maps, each spanned by Right
Ascension or Galactic longitude and Declination or Galactic latitude. Each
map corresponds to a specific photon energy. The energy binning of the map
cube may be either linear, logarithmic, or custom defined using an input
file.

:ref:`ctmapcube` generates a map cube FITS file comprising two extensions. The
primary extension contains a 3-dimensional image containing the map cube
values. A second extension named ``ENERGIES`` contains a binary table that
defines the energy values of each map in the cube.


General parameters
------------------

``inmodel [file]``
    Input model XML file.

``outcube [file]``
    Output map cube file.

``(ptsrcsig = 1.0) [real]``
    Sigma of Gaussian to be used for point sources (in arcmin).

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

``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.

``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.

``binsz [real]``
    Image scale (in degrees/pixel).

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
    Specifies whether the map cube should be published on VO Hub.

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

``(logfile = ctmapcube.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:ref:`csmodelsois`
