.. _csresspec:

csresspec
=========

Generates residual spectrum.


Synopsis
--------

This scripts generates a residual spectrum for a given model. It works for event
lists, counts cubes, On/Off observations or observation definition files. For
event lists parameters that define the spectral binning need to be provided so
that the script can bin the data internally. The model is then convolved with
the instrumental response function for the chosen binning (either specified by
the user or intrinsic to the observations) and used for residual computation. It
is possible to compute the model contributions from individual model components.

If starting from event lists or count cubes, before residual computation the 3D
counts and model cubes are collapsed into spectra by summing over all spatial
direction. It is possible to specify a mask to compute counts, model, and
residuals only for a region of interest.

For On/Off observations, counts, model, and residual spectra for the Off regions
are computed as well.

For an observation definition file, when all observations are event lists or
On/Off observations it is possible to stack them before residual computation. To
do so spatial binning parameters for event lists are required. Otherwise the
residuals are calculated separately for each observation.

Four options exist then for residual computation:

- the subtraction of the model from the counts (``SUB``)
- the subtraction and division by the model (``SUBDIV``)
- the subtraction and division by the square root of the model (``SUBDIVSQRT``)
- the computation of the significance based on the likelihood-ratio
  test for Poisson statistics (``SIGNIFICANCE``)

The ``SUBDIVSQRT`` and ``SIGNIFICANCE`` algorithms approximate the residual
significance and may become inaccurate in the low-counting regime.

The counts, model, and residual spectra are written into a FITS file.


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

``outfile [file]``
    Output residual spectrum file.

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.

``(components = no) [boolean]``
    Calculate model for individual components?

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

``stack [boolean]``
    Stack multiple observations for residual calculation?

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic) for stacking
    multiple event lists.

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method for stacking multiple event lists.

``xref [real]``
    Right Ascension / Galactic longitude of cube centre (J2000, in degrees)
    for stacking multiple event lists.

``yref [real]``
    Declination / Galactic latitude of cube centre (J2000, in degrees) for
    stacking multiple event lists.

``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude for stacking
    multiple event lists.

``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude for stacking
    multiple event lists.

``binsz [real]``
    Cube bin size (in degrees/pixel) for stacking multiple event lists.

``mask [boolean]``
    Mask data to calculate residuals in ROI?

``ra [real]``
    Right Ascension of circular selection region centre (J2000, in degrees).

``dec [real]``
    Declination of circular selection region centre (J2000, in degrees).

``rad [real]``
    Radius of circular selection region (in degrees).

``regfile [file]``
    Input exclusion region file in ds9 format.

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

``(logfile = csresspec.log) [string]``
    Log filename.


Related tools or scripts
------------------------

None
