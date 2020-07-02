.. _csresmap:

csresmap
========

Generates residual map.


Synopsis
--------

This scripts generates a residual map for a given model. It works for
event lists, counts cubes or observation definition files. For event
lists, parameters that define the spatial and spectral binning need to
be provided so that the script can bin the data internally. The model
is then convolved with the instrumental response function for that
binning and used for residual computation. Before residual computation,
the counts and model cubes are collapsed into maps by summing over all
energies. Four options exist then for residual computation:

- the subtraction of the model from the counts (``SUB``)
- the subtraction and division by the model (``SUBDIV``)
- the subtraction and division by the square root of the model (``SUBDIVSQRT``)
- the computation of the significance based on the likelihood-ratio
  test for Poisson statistics (``SIGNIFICANCE``)

The ``SUBDIVSQRT`` and ``SIGNIFICANCE`` algorithms approximate the residual
significance and may become inaccurate in the low-counting regime.

The residual map is written into a FITS file.  


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``modcube [file]``
    Input model cube file (generated with ctmodel).

``expcube [file]``
    Input exposure cube file.

``psfcube [file]``
    Input PSF cube file.

``edispcube [file]``
    Input energy dispersion cube file.

``bkgcube [file]``
    Input background cube file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``(edisp = no) [boolean]``
    Apply energy dispersion?

``outmap [file]``
    Output residual counts map file.

``(ebinalg = LOG) <FILE|LIN|LOG|POW> [string]``
    Algorithm for defining energy bins. For ``FILE``, the energy bins are defined
    in a FITS file that is specified by the ``ebinfile`` parameter, for ``LIN``
    ``LOG`` and ``POW`` there will be ``enumbins`` energy bins spaced linearly,
    logarithmically, or following a power law between ``emin`` and ``emax``,
    respectively. For ``POW``, the parameter ``ebingamma`` specifies the slope
    of the power law.

``emin [real]``
    Lower energy value for residual map (in TeV).

``emax [real]``
    Upper energy value for residual map (in TeV).

``(enumbins = 20) [integer]``
    Number of model cube energy bins for internal residual map computation.

``(ebinfile = NONE) [file]``
    Name of the file containing the energy binning definition if ``ebinalg=FILE``.
    You may use :ref:`csebins` to generate a file with appropriate energy binning.

``(ebingamma = 1.0) [real]``
    Exponent of the power law for ``POW`` energy binning. An exponent of 1.0
    corresponds to a logarithmic energy binning.

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of cube centre (J2000, in degrees).

``yref [real]``
    Declination / Galactic latitude of cube centre (J2000, in degrees).

``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.

``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.

``binsz [real]``
    Cube bin size (in degrees/pixel).

``algorithm <SUB|SUBDIV|SUBDIVSQRT|SIGNIFICANCE> [string]``
    Algorithm used to generate the residual map:

     ``SUB``: :math:`DATA - MODEL`

     ``SUBDIV``: :math:`(DATA - MODEL) / MODEL`

     ``SUBDIVSQRT``: :math:`(DATA - MODEL) / \sqrt{MODEL}`

     ``SIGNIFICANCE``: :math:`{\rm sign}(DATA-MODEL) \times \sqrt{ 2
     \times ( DATA \times \ln \left(\frac{DATA}{MODEL} \right) +
     MODEL - DATA ) }`


Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the residual map should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing residual map file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = csresmap.log) [string]``
    Log filename.


Related tools or scripts
------------------------

None
