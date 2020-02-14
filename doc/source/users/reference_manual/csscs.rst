.. _csscs:

csscs
==============================

Performs spectral component separation.


Synopsis
--------

This script performs a spectral component separation, that is,
produces a flux skymap of one or multiple sources by running a likelihood
analysis over multiple adjacent circular regions of interest (ROI)
based on prior knowledge of the sources' spectra. The centre of the
ROIs is displaced to coincide with the centre of each pixel in the
output flux map, while the radius is set through the user parameter ``rad``. 

The model fit per ROI is performed using :doc:`ctlike`, and only the
normalisations of the model components are treated as free parameters.
The script provides the possibility to fix sources other than the
source(s) of interest (hidden parameter ``fix_srcs``) or to fix the
normalisation of the background
model component(s) (hidden parameter ``fix_bkg``).

The script computes the source flux and its uncertainty in each
ROI/map pixel, as well as the significance of the source detection
(for faster execution the significance computation can be
deactivated using the hidden parameter ``calc_ts``). Optionally, the
script also computes an upper flux limit that is particularly useful
in case that the source is not significantly detected within an
ROI (hidden parameter ``calc_ulim``).

The script accepts in input event lists, binned data, or observation
definition XML files. Only data
within an energy interval spanned by ``emin`` and ``emax`` are
considered. For event lists two analysis methods are proposed:
unbinned or On/Off.

For the On/Off method  :doc:`csphagen` is
used to prepare On/Off observations for every ROI
analysed. In On/Off mode the user is required to provide an exclusion
map to avoid using regions of significant gamma-ray emission in the
evaluation of the background (exclusion maps can be generate using
:doc:`ctskymap`). Furthemore, On/Off analysis including
sources other than the source(s) of interest is not supported,
and if more than one source of interest is
present the response will be calculated based on the
hypothesis that emission from the sources of interest is uniform over
each ROI.

:ref:`csscs` generates a FITS file with empty primary extension.
Additional extensions provide for every source of interest a skymap of
flux, flux error, and, if requested, detection signficance and flux
upper limits.

General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcnames [string]``
    Semicolon-separated list of names of the sources in the source
    model XML file which should be considered for the spectral 
    component separation.

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

``inexclusion [file]``
    Optional FITS file containing a WCS map in the first extension that defines
    sky regions not to be used for background estimation (where map value != 0).

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation?

``outfile [file]``
    Output file.

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
    Declination / Galactic latitude of image centre (J2000, in
    degrees).

 ``rad [real]``
    Radius of region of interest for component separation (deg). Sets
    the correlation scale between neighbour pixels in the output maps. Must
    be at least sqrt(2) times ``binsz`` for full coverage of input data.

 ``emin [real]``
    Minimum energy (in TeV).

``emax [real]``
    Maximum energy (in TeV).

``method  <UNBINNED|ONOFF> [string]``
    For input event lists selects between ``UNBINNED`` analysis
    (3D spatial/energy likelihood without binning) and ``ONOFF``
    analysis (1D likelihood with background from Off regions).

``enumbins [integer]``
    Number of energy bins per light curve bin (for On/Off analysis only).

``(bkgmethod = REFLECTED) [string]``
    Method for background estimation in On/Off analysis.
    ``REFLECTED:`` background evaluated in regions with the same shape as
    source region reflected w.r.t. pointing direction for each
    observation.

``(srcshape = CIRCLE) [string]``
    Shape of the source region for On/Off analysis.
    ``CIRCLE``: circular region around given position.

``(bkgregmin = 2) [integer]``
    Minimum number of background regions that are required for an observation in
    On/Off analysis. If this number of background regions is not available the
    observation is skipped.

 ``(bkgregskip = 1) [integer]``
    Number of background regions that should be skipped next to the On regions.
    Typically, one region is skipped so that the Off regions are taken sufficiently
    distant from the On region, but in some cases it may be useful to keep the
    background regions next to the On region.

``(use_model_bkg = yes) [boolean]``
    Specifies whether the background model should be used for the computation
    of the ``alpha`` parameter and the predicted background rate in the Off
    region that is stored in the ``BACKRESP`` column of the Off spectrum when
    using the ``ONOFF`` method.

    If the parameter is set to ``no`` the background model is not used and the
    background rate is assumed identical within the On and Off regions. This
    is the classical IACT analysis method that is used when using reflected Off
    regions. In that case the ``alpha`` parameter becomes independent of energy
    and only reflects the ratio between the solid angles of the On and Off
    regions. The ``BACKRESP`` column in the Off spectrum will be filled with
    the solid angle of the On region. The data need to be fitted with the ``wstat``
    statistic, fitting with ``cstat`` will not work.

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
    all observations, ``WSTAT`` uses the W statistic for On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for
    binned or stacked observations.

``(calc_ts = yes) [boolean]``
    Compute TS value for each map bin?

``(calc_ulim = yes) [boolean]``
    Compute upper limit for each map bin?

``(fix_srcs = yes) [boolean]``
    Fix other sky model parameters?

``(fix_bkg = no) [boolean]``
    Fix background model parameters?


Standard parameters
-------------------

``(nthreads = 0) [integer]``
    Number of parallel processes (0=use all available CPUs).

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = csscs.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctlike`
:doc:`ctskymap`
:doc:`csphagen`
