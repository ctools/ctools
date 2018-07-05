.. _csspec:

csspec
======

Computes spectrum for a given source.


Synopsis
--------

This script computes the source spectrum by fitting a model in a given set
of spectral bins. The model fit per spectral bin is performed using :doc:`ctlike`
and the script provides the possibility to fix sources other than the
source of interest (hidden parameter ``fix_srcs``) or to fix the background
model component(s) (hidden parameter ``fix_bkg``). The script computes the
source flux and its uncertainty in each spectral bin, as well as the
significance of the source detection. Optionally, it also computes an upper
flux limit that is particularly useful in case that the source is not
significantly detected within a spectral bin (hidden parameter ``calc_ulim``).

There are two fundamental methods to run the script depending on the value of
the ``method`` parameter. If ``method=SLICE`` the energy interval defined by the
``emin`` and ``emax`` parameters is divided into a number of energy bins and an
independent maximum likelihood fit is performed in each of the energy bins,
while if ``method=NODES`` the spectral model will be replaced by a node function
that will be fit to all data. The latter is particularly useful if non-CTA data
should be fitted such as data from Fermi/LAT or COMPTEL. By default, ``method`` is
set to ``AUTO`` which will automatically select the method based on the input data.
For CTA-only observations ``SLICE`` will be used, otherwise ``NODES`` will be used.

The spectral binning is either defined by a FITS file containing the energy
boundaries of each bin (option ``ebinalg=FILE``) or as ``enumbins`` bins spread
linearly  (option ``ebinalg=LIN``) or logarithmically (option ``ebinalg=LOG``)
from a minimum energy, given by ``emin``, to a maximum energy, given by ``emax``.

For binned, stacked or On/Off CTA data, all energy bins that overlap with the
energy range spanned by ``emin`` and ``emax`` are considered. The number of spectral
bins is only approximately determined by the ``enumbins`` parameter. Naturally,
:ref:`csspec` cannot create more spectral bins than the number of energy bins that
are available in the data. In case that there are more energy bins in the data
than the number of spectral bins that are requested, :ref:`csspec` will fit the
data in all energy bins that overlap with a given spectral bin simultaneously.

:ref:`csspec` supports multiprocessing for ``method=SLICE``. By default the
analysis in each energy bin will be performed in parallel over as many processes
as the number of CPUs available on your machine. The maximum number of parallel
processes can be set by the user through the ``nthreads`` hidden parameter.

On output, the script will provide a FITS file with the fitted source 
spectrum in form of a binary table. Each row corresponds to a spectral bin.
The columns are the mean as well as the boundaries of the spectral bin, 
the fitted flux and flux error, the Test Statistics value (option
``calc_ts=yes``), the upper flux limit (option ``calc_ulim=yes``) and the
predicted number of events (only for unbinned data).

.. warning::
   The upper limit computation is not yet implemented for the ``NODES`` method.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for spectrum generation.

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
    Instrumental response function.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation?

``outfile [file]``
    Output spectrum FITS file.

``method <SLICE|NODES|AUTO> [string]``
    Spectrum generation method.
    ``SLICE`` will slice energy interval into ``enumbins`` independent bins,
    ``NODES`` will replace spectral model by node function with ``enumbins``
    nodes, and ``AUTO`` will automatically select the method based on the input
    data. For CTA-only observations ``SLICE`` will be used, otherwise ``NODES``
    will be used.

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

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.

``(calc_ts = yes) [boolean]``
    Compute TS for each spectral point?

``(calc_ulim = yes) [boolean]``
    Compute upper limit for each spectral point?

``(fix_srcs = yes) [boolean]``
    Fix other sky model parameters?

``(fix_bkg = no) [boolean]``
    Fix background model parameters?


Standard parameters
-------------------

``(nthreads = 0) [integer]``
    Number of parallel processes (0=use all available CPUs).

``(publish = no) [boolean]``
    Specifies whether the spectrum should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing source spectrum output file should be
    overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = csspec.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
