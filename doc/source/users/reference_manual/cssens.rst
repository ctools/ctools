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
using :doc:`ctlike` to determine it's detection significance as a function of
source flux. The source flux is then varied until the source significance
achieves a given level, specified by the (hidden) significance parameter
``sigma``.

The significance is estimated using the Test Statistic value, defined as twice
the log-likelihood difference that is obtained when fitting the simulated data
with and without the test source. The simplified assumption is made that the
significance (in Gaussian sigma) is the square root of the Test Statistic.

:ref:`cssens` will generate an ASCII file in comma-separated value (CSV) format
containing the sensitivity as function of energy. The first row is a header row
providing the column names. The following rows provide the mean logarithmic
energy and the boundaries of the energy bin for which the sensitivity was
computed. They also provide the flux within the energy bin in Crab units, in
photons (ph/cm2/s) and in energy (erg/cm2/s). Finally, the sensitivity is given
as the test source spectrum evaluated at the mean logarithmic energy multiplied
by the energy squared (erg/cm2/s).

The sensitivity can be displayed using the ``show_sensitivity.py`` script in the
example folder. Matplotlib is required to execute the script.


General parameters
------------------

``(inobs = NONE) [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for sensitivity computation.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation?

``(deadc = 0.98) [real]``
    Average deadtime correction factor.

``outfile [file]``
    ASCII file containing the sensitivity values.

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
    Number of energy bins (0=unbinned).

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
