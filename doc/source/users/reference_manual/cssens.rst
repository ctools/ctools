.. _cssens:

cssens
======

Computes differential or integrated CTA sensitivity.


Synopsis
--------

This script computes the differential or integrated CTA sensitivity using
maximum likelihood fitting of a test source. The differential sensitivity is
determined for a number of energy bins, the integral sensitivity is determined
for a number of energy thresholds. The test source is fitted to simulated data
using :doc:`ctlike` to determine its detection significance as a function of
source flux. The source flux is then varied until the source significance
achieves a given level, specified by the (hidden) significance parameter
``sigma``. As test source, any source in the input model definition XML file
specified by ``inmodel`` can be used; the test source name is specified by the
``srcname`` parameter.

The detection significance is estimated using the Test Statistic value, defined
as twice the log-likelihood difference that is obtained when fitting the
simulated data with and without the test source. The simplified assumption is
made that the significance (in Gaussian sigma) is the square root of the Test
Statistic.

By default the script simulates a single pointing with an exposure time
specified by the ``duration`` parameter. The pointing direction will be set to
the position of the test source, offset by the value specified by the ``offset``
parameter. If the test source has no position (which is for example the case
for a test source comprised of a map or map cube), a pointing direction of
(0,0) will be assumed in Right Ascension and Declination.

Alternatively, an observation definition XML file can be specified using the
hidden ``inobs`` parameter. In that way the user has full control over the
pointing sequence that should be simulated for the sensitivity estimation.

:ref:`cssens` supports multiprocessing. By default the analysis for each energy
bin/threshold will be performed in parallel over as many processes as the number of
CPUs available on your machine. The maximum number of parallel processes can be set
by the user through the ``nthreads`` hidden parameter.

:ref:`cssens` will generate a FITS file with a single extension that contains a
binary table with the sensitivity as function of energy.

The sensitivity FITS file can be displayed using the ``show_sensitivity.py`` script
in the example folder. Matplotlib is required to execute the script.


General parameters
------------------

``(inobs = NONE) [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for sensitivity computation.

``(instrument = CTA) [string]``
    Name of Cherenkov telescope.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation?

``(deadc = 0.98) [real]``
    Average deadtime correction factor.

``outfile [file]``
    Output sensitivity FITS file.

``(offset = 0.0) [real]``
    Offset angle of source in field of view (in degrees).	 

``duration [real]``
    Effective exposure time (in seconds).

``rad [real]``
    Radius of Region of Interest (RoI) (in degrees).

``emin [real]``
    Lower energy limit for sensitivity computation (in TeV).

``emax [real]``
    Upper energy limit for sensitivity computation (in TeV).

``bins [integer]``
    Number of energy bins for sensitivity computation.

``(enumbins = 0) [integer]``
    Number of energy bins for binned analysis (``0`` = unbinned analysis).

``(npix = 200) [integer]``
    Number of pixels for binned analysis.

``(binsz = 0.05) [real]``
    Pixel size for binned analysis.

``(type = Differential) <Differential|Integral> [string]``
    Sensitivity type.

``(sigma = 5.0) [real]``
    Significance threshold.

``(max_iter = 50) [integer]``
    Maximum number of iterations.

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.

``(mincounts = 10) [integer]``
    Constraint on the minimum number of required source counts. Conventionally,
    a constraint for a minimum number of 10 source counts is applied for CTA
    sensitivity estimates. If ``0`` is specified then no source counts limit
    will be applied.

``(bkgexcess = 0.0) [real]``
    Constraint on the minimum number of required source counts with respect to
    the number of background counts. This value is a fraction that is
    conventionally set to ``0.05`` (or 5%) which means that the number of source
    counts needs to exceed 5% of the number of background counts. If ``0.0`` is
    specified then no constraint will be applied.

``(bkgrad = 0.33) [real]``
    Radius in degrees used to estimate the number of background counts
    underlying the source. This radius is only used if ``bkgexcess > 0.0``.


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
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = cssens.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
